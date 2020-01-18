"""Testing script for offline analysis.

We consider mixed atoms, and compare the paths taken between the observations.

We iterate through function/property pairs, then group pairs of observations by the pairs of instrumentation points."""
import sys
import pprint

# let Python find VyPR
sys.path.append("../local_test")

import VyPRAnalysis
from VyPR.SCFG.parse_tree import ParseTree

# set configuration with the default file at config.json
VyPRAnalysis.set_config_file()

if __name__ == "__main__":
    functions = VyPRAnalysis.list_functions()
    for function in functions:
        print("Processing function %s" % function)
        # construct the SCFG to be used in path comparison
        function_scfg = function.get_scfg()
        # derive the context-free grammar to be used in path comparison
        grammar = function_scfg.derive_grammar()
        VyPRAnalysis.utils.write_scfg(function_scfg, 'gvs/%s-scfg.gv' % function.fully_qualified_name)
        calls = function.get_calls()
        # construct a map from pairs of instrumentation points to pairs of observations
        instrumentation_point_pair_map = {}
        # for each call of this function, we iterate through the verdicts and find the collapsing observations
        for call in calls:
            # get the verdicts during this call for the property
            verdicts = call.get_verdicts()
            for verdict in verdicts:
                collapsing_atom_index = verdict.collapsing_atom
                # check that the atom that caused collapse is a normal atom, and not mixed
                atom = VyPRAnalysis.Atom(property_hash=function.property, index_in_atoms=collapsing_atom_index)
                if VyPRAnalysis.utils.get_atom_category(atom.get_structure()) == "mixed":
                    # get the observation for this verdict/index combination
                    observations = list(
                        filter(
                            lambda obs: obs.atom_index == collapsing_atom_index, verdict.get_observations()
                        )
                    )
                    # sort by atom sub index
                    observations = sorted(observations, key=lambda observation : observation.observation_time)

                    # compute the path between the observations
                    early_path = observations[0].reconstruct_reaching_path(function_scfg)
                    late_path = observations[1].reconstruct_reaching_path(function_scfg)
                    # compute the difference
                    difference_path = late_path[len(early_path):]

                    # we will have 2 observations, so extract the pair of instrumentation points
                    # our goal is to construct a map from instrumentation point pairs to observation pairs
                    instrumentation_point_pair = (observations[0].instrumentation_point,
                                                  observations[1].instrumentation_point)

                    # compute the path between the two observations

                    if instrumentation_point_pair_map.get(instrumentation_point_pair):
                        if instrumentation_point_pair_map[instrumentation_point_pair].get(verdict.verdict):
                            instrumentation_point_pair_map[instrumentation_point_pair][verdict.verdict]\
                                .append(difference_path)
                        else:
                            instrumentation_point_pair_map[instrumentation_point_pair][verdict.verdict]\
                                = [difference_path]
                    else:
                        instrumentation_point_pair_map[instrumentation_point_pair]\
                            = {verdict.verdict : [difference_path]}

        # iterate through the map from instrumentation point pairs to verdicts to difference paths
        # and compute the intersections of the difference paths for each verdict
        for instrumentation_point_pair in instrumentation_point_pair_map:
            for verdict in instrumentation_point_pair_map[instrumentation_point_pair]:
                # replace the list of difference paths with a single parametric path
                # for this, we need to compute the intersection parse tree of all the parse trees,
                # and then read off its leaves
                paths = instrumentation_point_pair_map[instrumentation_point_pair][verdict]
                first_parse_tree = ParseTree(paths[0], grammar, paths[0][0]._source_state)
                other_parse_trees = list(
                    map(
                        lambda path: ParseTree(path, grammar, path[0]._source_state),
                        paths[1:]
                    )
                )
                intersection = first_parse_tree.intersect(other_parse_trees)
                intersection.write_to_file("gvs/inst-pair-%i-%i-intersection.gv" %
                                           instrumentation_point_pair)
                parametric_path = intersection.read_leaves()
                # replace the list of paths with the intersection
                instrumentation_point_pair_map[instrumentation_point_pair][verdict] = parametric_path

        pprint.pprint(instrumentation_point_pair_map)