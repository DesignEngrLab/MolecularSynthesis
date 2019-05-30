using GraphSynth.Representation;
using System;
using System.IO;
using OpenBabelFunctions;
using GraphSynth.Search.Tools;



namespace GraphSynth.Search.Algorithms {

    /// <summary>
    /// Randomly select a valid option at any step.
    /// </summary>
    public class Random : AbstractAlgorithm
    {

        public Random(GlobalSettings settings) : base(settings)
        {
            
        }


        public option ChooseOption(candidate cand)
        {
            var options = GetAvailableOptions(cand);
            return options.Count > 0 ? options[Rand.Next(options.Count)] : null;
        }

        public option ChooseCarboxOption(candidate cand)
        {
            var options = GetCarboxylOptions(cand);
            return options.Count > 0 ? options[Rand.Next(options.Count)] : null;
        }

        public option ChooseCarboxOptionBestAngle(candidate cand)
        {
            var options = GetCarboxylOptions(cand);
            option bestOpt = null;
            var bestAngle = .0;
            foreach (var opt in options) 
            {
                var evalcand = CopyAndApplyOption(opt, cand, true);
                var mol = OBFunctions.designgraphtomol(evalcand.graph);
                var angle = CalAngle(mol);
                if (angle > 180)
                {
                    Console.WriteLine(angle + " too large");
                    Environment.Exit(0);
                }
                if (angle > bestAngle)
                {
                    bestAngle = angle;
                    bestOpt = opt;
                }
            }
            return bestOpt;
        }

        public option ChooseCarboxOptionUsingEstimator(candidate cand, Computation cpt, MessageClient clt)
        {
            var options = GetCarboxylOptions(cand);
            option bestOpt = null;
            var bestProperty = .0;
            foreach (var opt in options) 
            {
                var evalcand = CopyAndApplyOption(opt, cand, true);
                var mol = OBFunctions.designgraphtomol(evalcand.graph);
                var linkerName = GetLinkerName(evalcand);
                var coeff = Path.Combine("possible", "linker" + linkerName + ".coeff");
                var lmpdat = Path.Combine("possible", "linker" + linkerName + ".lmpdat");
                Converter.moltoUFF(OBFunctions.designgraphtomol(cand.graph), coeff, lmpdat, false, 100);
                computation.CalculateFeature(linkerName, true);
                Environment.Exit(0);

                if (property > bestProperty)
                {
                    bestProperty = property;
                    bestOpt = opt;
                }
            }
            
        }


    }



}