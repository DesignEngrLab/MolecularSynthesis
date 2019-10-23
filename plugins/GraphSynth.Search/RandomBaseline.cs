using System;
using System.Collections.Generic;
using System.IO;
using GraphSynth.Representation;
using GraphSynth.Search.Algorithms;
using GraphSynth.Search.Tools;
using OpenBabelFunctions;
namespace GraphSynth.Search
{
    class RandomBaseline : SearchProcess
    {
        private readonly string _runDirectory;

        private candidate Seed;
        private JobBuffer jobBuffer;
        private Computation computation;
        private StreamWriter writer;

        private const double PROB_DEC = 0.1;
        private const double PROB_NON_TEM_INIT = 1;
        private const double PROB_TEM_INIT = 1 - PROB_NON_TEM_INIT;
        private const int NUM_EPOCH = 20;
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
            Console.WriteLine("Fall 2019 One Drive....");
            for (var e = 0; e < NUM_EPOCH; e++)
            {
                var rand = rnd.NextDouble();
                Console.WriteLine("{}", rand);
            }


        }

        //private candidate step(candidate cand)
        //{
        //    var randFloat = 
        //    return cand;
        //}



        public override string text => "RandomBaseline Search Runner";

    }


}
