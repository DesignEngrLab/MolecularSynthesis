using GraphSynth.Representation;
using System;
using OpenBabelFunctions;



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
            foreach (var opt in options) 
            {
                var evalcand = CopyAndApplyOption(opt, cand, true);
                var mol = OBFunctions.designgraphtomol(evalcand.graph);
                var angle = CalAngle(mol);
                Console.WriteLine(angle);
            }
            Environment.Exit(0);
            return options.Count > 0 ? options[Rand.Next(options.Count)] : null;
        }


    }



}