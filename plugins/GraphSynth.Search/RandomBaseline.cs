using System;
using System.IO;
using GraphSynth.Representation;
using GraphSynth.Search.Algorithms;
using GraphSynth.Search.Tools;
using OpenBabelFunctions;


namespace GraphSynth.Search
{
    public class RandomBaseline : SearchProcess
    {
        private readonly string _runDirectory;

        private candidate Seed;
        private JobBuffer jobBuffer;
        private Computation computation;
        private StreamWriter writer;
        
        private const int NUM_EPOCH = 100;
        private const int NUM_RUNS = 20;
        private System.Random rnd = new System.Random();
        
        // Baseline: With probability distribution [1-0.1*steps, 0.1*steps]
        // Choose to a non-terminal rule or a terminal rule
        // If any one is infeasible at current moment, choose the other. If both are infeasible, restart.
        public RandomBaseline(GlobalSettings settings) : base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 1;
            AutoPlay = true;

            Seed = new candidate(OBFunctions.tagconvexhullpoints(settings.seed), settings.numOfRuleSets);

            _runDirectory = Path.Combine(settings.OutputDirAbs, "RandomBaseline");
            if (Directory.Exists(_runDirectory))
                Directory.Delete(_runDirectory, true);
            Directory.CreateDirectory(_runDirectory);
            var learnDirectory = Path.Combine(settings.OutputDirAbs, "morfLearn");
            computation = new Computation(_runDirectory, learnDirectory, "point", "stiff");
            writer = new StreamWriter(Path.Combine(_runDirectory, "RandomBaseline.txt"));
            jobBuffer = new JobBuffer(_runDirectory);
        }

        protected override void Run()
        {
            Console.WriteLine("Fall 2019 One Drive New Mac Matt....");
            var agent = new Algorithms.Random(settings);

            for (var r = 0; r < NUM_RUNS; r++)
            {
                for (var e = 0; e < NUM_EPOCH; e++)
                {
                    Console.WriteLine("Epoch: {0}", e);
                    candidate cand = null;
                    while (cand == null)
                    {
                        cand = Seed.copy();
                        while (true)
                        {
                            cand = agent.ChooseAndApplyAnyOption(cand);
                            if (cand == null)
                            {
                                Console.WriteLine("Fail, Rebuild");
                                break;
                            }
                            if (agent.IsTerminalCandidate(cand))
                            {
                                Console.WriteLine("Choose terminal rule.");
                                break;
                            }
                            Console.WriteLine("Choose non-terminal rule.");
                        }
                    }
                    var linkerName = AbstractAlgorithm.GetLinkerName(cand);
                    Console.WriteLine(linkerName);
                    var atomCount = AbstractAlgorithm.CountAtoms(cand);
                    Console.WriteLine(atomCount);
                    var smi = OBFunctions.moltoSMILES(OBFunctions.designgraphtomol(cand.graph));
                    Console.WriteLine(smi);
                    writer.Write("{0} ", atomCount);
                }
                writer.WriteLine();
            }
            writer.Close();
        }
        
        public override string text => "RandomBaseline Search Runner";

    }
}
