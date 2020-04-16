using GraphSynth;
using GraphSynth.Representation;
using MolecularSynthesis.Tools;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;



namespace MolecularSynthesis.Algorithms {

    /// <summary>
    /// Randomly select a valid option at any step.
    /// </summary>
    public class RandomAlgorithm : AbstractAlgorithm
    {
        public RandomAlgorithm(GlobalSettings settings) : base(settings)
        {
             
        }
        
        
        public candidate ChooseAndApplyAnyOption(candidate cand)
        {
            var options = GetAllOptions(cand);
            return options.Count > 0 ? CopyAndApplyOption(options[Rand.Next(options.Count)], cand, true) : null;
        }
        
        public bool IsTerminalCandidate(candidate cand)
        {
            var numRule = cand.recipe.Count;
            return cand.recipe[numRule - 1].ruleSetIndex == 1;
        }
        

        public candidate ChooseAndApplyCarboxOptionUsingEstimator(candidate cand, Computation cpt, MessageClient clt, string runDir, int epoch)
        {
            option bestOpt = null;
            var bestProperty = Double.NegativeInfinity;
            var options = GetTerminalOptions(cand);
            foreach (var opt in options) 
            {
                var evalcand = CopyAndApplyOption(opt, cand, true);
                var mol = OBFunctions.designgraphtomol(evalcand.graph);
                var linkerName = AbstractAlgorithm.GetLinkerName(evalcand) + "-E" + epoch.ToString();
                var coeff = Path.Combine(runDir, "data", "linker" + linkerName + ".coeff");
                var lmpdat = Path.Combine(runDir, "data", "linker" + linkerName + ".lmpdat");
                cpt.CalculateFeature(linkerName);
                var properpty = Convert.ToDouble(clt.SendMessage("[Predict]" + " " + linkerName));

                if (properpty > bestProperty)
                {
                    bestProperty = properpty;
                    bestOpt = opt;
                }
            }
            Debug.WriteLine("Best {0}", bestProperty);
            return bestOpt == null ? null : CopyAndApplyOption(bestOpt, cand, true);
        }


    }



}