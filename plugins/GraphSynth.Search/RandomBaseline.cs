using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using GraphSynth.Representation;
using GraphSynth.Search.Algorithms;
using GraphSynth.Search.Tools;
using OpenBabelFunctions;
using Random = System.Random;


namespace GraphSynth.Search
{
    public class RandomBaseline : SearchProcess
    {
        private readonly string _runDirectory;

        private candidate Seed;
        private JobBuffer jobBuffer;
        private Computation computation;
        private StreamWriter writer;
        private Algorithms.Random agent;

        private const int GEN_SET_SIZE = 10;
        private const int NUM_EPOCH = 100;
        private const int NUM_RUNS = 20;

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
            agent = new Algorithms.Random(settings);
        }

        protected override void Run()
        {
            
            Console.WriteLine("Fall 2019 One Drive New Mac Matt....");
            for (var r = 0; r < NUM_RUNS; r++)
            {
                Dictionary<string, int> MolSet = new Dictionary<string, int>();
                for (var s = 0; s < GEN_SET_SIZE; s++)
                {
                    var cand = GenerateCand();
                    var atomNum = AbstractAlgorithm.CountAtoms(cand);
                    var smi = OBFunctions.moltoSMILES(OBFunctions.designgraphtomol(cand.graph));
                    Console.WriteLine(smi);
                    if (!MolSet.ContainsKey(smi))
                        MolSet.Add(smi, atomNum);
                }
                for (var e = 0; e < NUM_EPOCH; e++)
                {
                    Console.WriteLine("Epoch: {0}", e);
                    var cand = GenerateCand();
                    var atomNum = AbstractAlgorithm.CountAtoms(cand);
                    var smi = OBFunctions.moltoSMILES(OBFunctions.designgraphtomol(cand.graph));
                    Console.WriteLine(smi);
                    if (!MolSet.ContainsKey(smi))
                        MolSet.Add(smi, atomNum);
                    writer.Write("{0} ", MolSet.Values.ToList().Max());
                }
                writer.WriteLine();
            }
            writer.Close();
        }

        private candidate GenerateCand()
        {
            candidate cand = null;
            while (cand == null)
            {
                cand = Seed.copy();
                while (true)
                {
                    try
                    {
                        cand = agent.ChooseAndApplyAnyOption(cand);
                        if (cand == null)
                        {
                            //Console.WriteLine("Fail, Rebuild");
                            break;
                        }

                        if (agent.IsTerminalCandidate(cand))
                        {
                            //Console.WriteLine("Choose terminal rule.");
                            break;
                        }
                        //Console.WriteLine("Choose non-terminal rule.");
                    }
                    catch (Exception e)
                    {
                        Console.WriteLine("Convex Hull Fails.");
                        cand = Seed.copy();
                    }
                }
            }
            return cand;
        }

    public override string text => "RandomBaseline Search Runner";

    }
}
