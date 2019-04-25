using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using GraphSynth.Representation;
using GraphSynth.Search.Algorithms;
using GraphSynth.Search.Tools;
using OpenBabelFunctions;





namespace GraphSynth.Search
{
    public class RandomRuleApplication: SearchProcess
    {
        
        private readonly string _runDirectory;
        private readonly string _inputFilePath;
        
        private candidate Seed;
        private const int NUM_TRAIL = 5;
        private const int TOTAL_RULE = 5;


        /// <inheritdoc />
        /// <summary>
        /// Initializes SearchProcess properties.
        /// </summary>
        public RandomRuleApplication(GlobalSettings settings): base(settings) 
        {
            RequireSeed = true;
            RequiredNumRuleSets = 1;
            AutoPlay = true;
            
            _runDirectory = Path.Combine(settings.OutputDirAbs, "RandomRuleApplication");
            if (!Directory.Exists(_runDirectory))
                Directory.CreateDirectory(_runDirectory);
            Seed = new candidate(OBFunctions.tagconvexhullpoints(settings.seed), settings.numOfRuleSets);
        }

        protected override void Run()
        {
            var linkerSet = new HashSet<string>();
            var linkerBuffer = new MolBuffer(_runDirectory);
            
            var agent = new Algorithms.Random(settings);
            for (var t = 0; t < NUM_TRAIL; t++)
            {
                Console.WriteLine("Trail: {0}", t);
                var cand = Seed.copy();

                for (var step = 1; step < TOTAL_RULE + 1; step++)
                {
                    var opt = agent.ChooseOption(cand);
                    if (opt == null)
                        return;
                    AbstractAlgorithm.ApplyOption(opt, cand, true);
                }
                var carboxOpt = agent.ChooseCarboxOption(cand);
                if (carboxOpt == null)
                    return;
                AbstractAlgorithm.ApplyOption(carboxOpt, cand, true);
                var candSMILE = OBFunctions.moltoSMILES(OBFunctions.designgraphtomol(cand.graph));
                var linkerName = AbstractAlgorithm.GetLinkerName(cand);
                Console.WriteLine(candSMILE);
                Console.WriteLine(linkerName);
                if (linkerSet.Contains(linkerName))
                {
                    t--;
                    continue;
                }
                linkerSet.Add(linkerName);
                linkerBuffer.Add(linkerName, cand, AbstractAlgorithm.Rand.NextDouble());
            }

            while (linkerBuffer.Size() > 0)
            {
                linkerBuffer.Remove();
            }

        }
        
        public override string text => "RandomTrail Search Runner";
        
    }
}