"""Testing script for offline analysis.

We iterate through function/property pairs, then group observations and the function calls occurring during them
(so below in the call chain) with respect to instrumentation points.

We then compare the paths taken by the function calls lower down in the call down for observations
that resulted in failure vs observations that resulted in success.

This gives us an indication of whether paths taken through functions that are called during observations
may be responsible for an observation violating a constraint.

We construct a map instrumentation point -> verdict -> function -> {call paths}, using which we can compare paths from
the same function for a fixed instrumentation point and verdict."""
import sys
import pprint

# let Python find VyPR
sys.path.append("../local_test")

import VyPRAnalysis

# set configuration with the default file at config.json
VyPRAnalysis.set_config_file()

import VyPRAnalysis.interprocedural
import VyPRAnalysis.utils

from VyPR.SCFG.parse_tree import ParseTree

if __name__ == "__main__":
    VyPRAnalysis.set_server("http://localhost:9002/")
    functions = VyPRAnalysis.list_functions()
    function_id_to_scfg = {}
    for function in functions:
        # construct this function's SCFG - we need this for path reconstruction
        # note, for path comparison we will need to use the same SCFG (wih the same address in memory) throughout
        if not (function.id in function_id_to_scfg):
            scfg = function.get_scfg()
            function_id_to_scfg[function.id] = scfg
        else:
            scfg = function_id_to_scfg[function.id]
        # VyPRAnalysis.utils.write_scfg(scfg, 'gvs/%s-scfg.gv' % function.fully_qualified_name)
        # we care about instrumentation points that generated observations that caused a verdict to be reached
        # we construct a map from instrumentation points IDs to a set of observations causing failure
        # and a set of observations causing success.
        instrumentation_point_id_map = {}
        calls = function.get_calls()
        for call in calls:
            transaction = VyPRAnalysis.transaction(id=call.trans)

            # get the verdicts during this call for the property
            verdicts = call.get_verdicts()

            for verdict in verdicts:
                collapsing_atom_index = verdict.collapsing_atom
                collapsing_atom_sub_index = verdict.collapsing_atom_sub_index
                # check that the atom that caused collapse is a normal atom, and not mixed
                atom = VyPRAnalysis.Atom(property_hash=function.property, index_in_atoms=collapsing_atom_index)
                if VyPRAnalysis.utils.get_atom_category(atom.get_structure()) == "normal":
                    # get the observation for this verdict/index combination
                    observation = list(
                        filter(
                            lambda obs: obs.atom_index == collapsing_atom_index and
                                        obs.sub_index == collapsing_atom_sub_index,
                            verdict.get_observations()
                        )
                    )[0]
                    inst_point_id = observation.instrumentation_point
                    if instrumentation_point_id_map.get(inst_point_id):
                        if instrumentation_point_id_map[inst_point_id].get(verdict.verdict):
                            instrumentation_point_id_map[inst_point_id][verdict.verdict].append(observation)
                        else:
                            instrumentation_point_id_map[inst_point_id][verdict.verdict] = [observation]
                    else:
                        instrumentation_point_id_map[inst_point_id] = {verdict.verdict: [observation]}

        # expand the structure we just created to map from observations to callees occurring during those observations
        for inst_point_id in instrumentation_point_id_map:
            for verdict_value in instrumentation_point_id_map[inst_point_id]:
                # we only care about the paths taken by a function, not with respect to observations
                function_id_map = {}
                for obs in instrumentation_point_id_map[inst_point_id][verdict_value]:
                    # get the transition of which this observation is a member
                    verdict = VyPRAnalysis.verdict(id=obs.verdict)
                    function_call = VyPRAnalysis.function_call(id=verdict.function_call)
                    transaction = VyPRAnalysis.transaction(id=function_call.trans)
                    # get the call tree of the transaction
                    call_tree = VyPRAnalysis.interprocedural.CallTree(transaction)
                    # get the callees of function_call
                    callees = call_tree.get_direct_callees(function_call)
                    # filter the callees to those that occurred during obs
                    callees = list(
                        filter(
                            lambda callee: obs.observation_time < callee.time_of_call < obs.observation_end_time,
                            callees
                        )
                    )
                    for callee in callees:
                        callee_function = VyPRAnalysis.function(id=callee.function)
                        if function_id_map.get(callee.function):
                            callee_function_scfg = function_id_to_scfg[callee.function]
                            path = callee.reconstruct_path(callee_function_scfg)
                            function_id_map[callee.function].append(path)
                        else:
                            if function_id_to_scfg.get(callee.function):
                                callee_function_scfg = function_id_to_scfg[callee.function]
                            else:
                                callee_function_scfg = callee_function.get_scfg()
                                function_id_to_scfg[callee.function] = callee_function_scfg
                            path = callee.reconstruct_path(callee_function_scfg)
                            function_id_map[callee.function] = [path]

                # reassign this part of instrumentation_point_id_map with a dictionary
                instrumentation_point_id_map[inst_point_id][verdict_value] = function_id_map

            # now, for this function/property pair, for each instrumentation point we compare paths taken
            # through each function that we found
            for inst_point_id in instrumentation_point_id_map:
                for verdict_value in instrumentation_point_id_map[inst_point_id]:
                    for function_id in instrumentation_point_id_map[inst_point_id][verdict_value]:
                        paths = instrumentation_point_id_map[inst_point_id][verdict_value][function_id]
                        relevant_scfg = function_id_to_scfg[function_id]
                        grammar = relevant_scfg.derive_grammar()
                        parse_tree = ParseTree(paths[0], grammar, relevant_scfg.starting_vertices)
                        # intersect with the other paths
                        other_parse_trees = list(
                            map(
                                lambda path: ParseTree(path, grammar, relevant_scfg.starting_vertices),
                                paths[1:]
                            )
                        )
                        intersection_parse_tree = parse_tree.intersect(other_parse_trees)
                        intersection_parse_tree.write_to_file("gvs/inst-%i-v-%i-f%i-intersection-parse-tree.gv" %
                                                              (function_id, verdict_value, function_id))

            print("For function/property pair %s" % function)
            pprint.pprint(instrumentation_point_id_map)