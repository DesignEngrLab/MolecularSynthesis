using System;
using System.Collections;
using System.IO;
using GraphSynth.Representation;
using GraphSynth.Search.Algorithms;
using OpenBabelFunctions;





namespace GraphSynth.Search
{
    public class RandomRuleApplication: SearchProcess
    {
        
        private string _runDirectory;
        private string _inputFilePath;
        
        private candidate Seed;
        private int NUM_TRAIL = 5;
        private int TOTAL_RULE = 5;


        /// <inheritdoc />
        /// <summary>
        /// Initializes SearchProcess properties.
        /// </summary>
        public RandomRuleApplication() {
            RequireSeed = true;
            RequiredNumRuleSets = 1;
            AutoPlay = true;
        }

        protected override void Run() {
            _runDirectory = Path.Combine(settings.OutputDirAbs, "RandomRuleApplication");
            if (!Directory.Exists(_runDirectory))
                Directory.CreateDirectory(_runDirectory);
            Seed = new candidate(OBFunctions.tagconvexhullpoints(settings.seed), settings.numOfRuleSets);
            
            
            var agent = new Algorithms.Random(settings);
            for (var t = 0; t < NUM_TRAIL; t++)
            {
                Console.WriteLine("Trail: {0}", t);
                var cand = Seed.copy();
                var linkerDir = _runDirectory + "/trail" + t + "/";
                Directory.CreateDirectory(linkerDir);
                var coeff = Path.Combine(linkerDir, "linker0.coeff");
                var lmpdat = Path.Combine(linkerDir, "linker0.lmpdat");
                // Set up UFF and run lammps
                AbstractAlgorithm.Converter.moltoUFF(OBFunctions.designgraphtomol(cand.graph), coeff, lmpdat, false, 100);

                var candSMILE = OBFunctions.moltoSMILES(OBFunctions.designgraphtomol(cand.graph));
                
                Console.WriteLine("{0}", candSMILE);
                for (var step = 1; step < TOTAL_RULE + 1; step++)
                {
                    var opt = agent.ChooseOption(cand);
                    if (opt == null)
                        return;
                    agent.ApplyOption(opt, cand, true);
                    var mol = OBFunctions.designgraphtomol(cand.graph);
                    candSMILE = OBFunctions.moltoSMILES(mol);
                    Console.WriteLine("{0}", candSMILE);
                    coeff = Path.Combine(linkerDir, "linker" + step + ".coeff");
                    lmpdat = Path.Combine(linkerDir, "linker" + step + ".lmpdat");

                    AbstractAlgorithm.Converter.moltoUFF(mol, coeff, lmpdat, false, 100);
                }
                var carboxOpt = agent.ChooseCarboxOption(cand);
                if (carboxOpt == null)
                    return;
                agent.ApplyOption(carboxOpt, cand, true);
                candSMILE = OBFunctions.moltoSMILES(OBFunctions.designgraphtomol(cand.graph));
                Console.WriteLine("{0}", candSMILE);
                coeff = Path.Combine(linkerDir, "linker" + (TOTAL_RULE+1) + ".coeff");
                lmpdat = Path.Combine(linkerDir, "linker" + (TOTAL_RULE+1) + ".lmpdat");

                AbstractAlgorithm.Converter.moltoUFF(OBFunctions.designgraphtomol(cand.graph), coeff, lmpdat, false, 100);
            }

        }
        
        public override string text => "RandomTrail Search Runner";
        
    }
}