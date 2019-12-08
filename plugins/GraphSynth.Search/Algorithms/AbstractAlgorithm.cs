using GraphSynth.Representation;
using OpenBabel;
using OpenBabelFunctions;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;



namespace GraphSynth.Search.Algorithms {
    
    public abstract class AbstractAlgorithm {
        public static GlobalSettings Settings;
        public static readonly System.Random Rand = new System.Random();
        private static readonly string IODir = OBFunctions.GetRamDir();
        private const double AngleFloor = 155; // minimum acceptable angle between carboxylates

        
        protected AbstractAlgorithm(GlobalSettings settings_) {
            Settings = settings_;
        }

        public static string GetLinkerName(candidate cand)
        {
            var arr = cand.recipe.Select(x => Convert.ToString(x.optionNumber));
            return String.Join("-", arr);
        }
        
        
        /// <summary>
        /// Get all available options for the given graph.
        /// </summary>
        /// <param name="cand"></param>
        /// <returns></returns>
        public static List<option> GetAllOptions(candidate cand) {
            var options0 = new List<option>();
            options0.AddRange(Settings.rulesets[0].recognize(cand.graph, false));
            var options1 = new List<option>();
            options1.AddRange(Settings.rulesets[1].recognize(cand.graph, false));
            return options0.Concat(options1).ToList();
        }
        

        /// <summary>
        /// Returns the option corresponding to the carboxyl rule application.
        /// </summary>
        /// <param name="cand"></param>
        /// <returns></returns>
        public static List<option> GetTerminalOptions(candidate cand) {
            var options = new List<option>();
            options.AddRange(Settings.rulesets[1].recognize(cand.graph, false));
            return options;
        }

        /// <summary>
        /// Apply the option to the candidate and store the agent's evaluation.
        /// </summary>
        public void ApplyOption(option opt, candidate cand, bool doMinimize) {
            cand.graph.globalVariables.Add(cand.f0); // track fitness values of previous states
            opt.apply(cand.graph, null);
            cand.addToRecipe(opt);
/*            if(doMinimize)
                cand.graph = Minimize(cand.graph);*/
            cand.graph = OBFunctions.tagconvexhullpoints(cand.graph);
        }


        /// <summary>
        /// Copies the candidate graph, transfers the L-mapping, and returns the resultant candidate.
        /// </summary>
        public candidate CopyAndApplyOption(option opt, candidate cand, bool doMinimize) {
            var newCand = cand.copy();
            var newOpt = opt.copy();
            SearchProcess.transferLmappingToChild(newCand.graph, cand.graph, newOpt);
            ApplyOption(newOpt, newCand, doMinimize);
            return newCand;
        }
        
        /// <summary>
        /// Clean way to minimize a graph.
        /// </summary>
        private designGraph Minimize(designGraph graph) {
            var mol = new OBMol();
            OBFunctions.updatepositions(graph, mol);
            return graph;
        }
        
    }
}