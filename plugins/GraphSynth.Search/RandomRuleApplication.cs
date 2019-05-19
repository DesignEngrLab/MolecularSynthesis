using System;
using System.Collections.Generic;
using System.IO;
using System.Threading;
using GraphSynth.Representation;
using GraphSynth.Search.Algorithms;
using GraphSynth.Search.Tools;
using OpenBabelFunctions;





namespace GraphSynth.Search
{
    public class RandomRuleApplication: SearchProcess
    {
        
        private readonly string _runDirectory;

        private candidate Seed;
        private JobBuffer jobBuffer;
        private LearningServer server;

        private const int NUM_TRAIL = 5;
        private const int TOTAL_RULE = 5;
        
        private static Mutex mutex = new Mutex();


        /// <inheritdoc />
        /// <summary>
        /// Initializes SearchProcess properties.
        /// </summary>
        public RandomRuleApplication(GlobalSettings settings): base(settings) 
        {
            RequireSeed = true;
            RequiredNumRuleSets = 1;
            AutoPlay = true;
            
            Seed = new candidate(OBFunctions.tagconvexhullpoints(settings.seed), settings.numOfRuleSets);
            
            _runDirectory = Path.Combine(settings.OutputDirAbs, "RandomRuleApplication", "randomCarbox");
            var learnDirectory = Path.Combine(settings.OutputDirAbs, "morfLearn");
            if (Directory.Exists(_runDirectory))
                Directory.Delete(_runDirectory, true);
            Directory.CreateDirectory(_runDirectory);

            jobBuffer = new JobBuffer(_runDirectory);
            server = new LearningServer(_runDirectory, learnDirectory, "point", "stiff");
        }

        protected override void Run()
        {
            Thread generateLinkers = new Thread(Generate);
            Thread autoReleaseBuffer = new Thread(AutoSubmitSimulation);
            
            generateLinkers.Start();
            autoReleaseBuffer.Start();

            generateLinkers.Join();
            autoReleaseBuffer.Join();

        }
        
        private void AutoSubmitSimulation()
        {
            var allSubmitted = false;
            while (true)
            {
                //mutex.WaitOne();
                if (jobBuffer.CanFeedIn())
                {
                    allSubmitted = jobBuffer.Simulate();
                }
                var allFinished = jobBuffer.Check_finised(server);
                if (allFinished && allSubmitted)
                    break;
                //mutex.ReleaseMutex();
            }
        }

        private void Generate()
        {
            Console.WriteLine(NUM_TRAIL);
            var linkerSet = new HashSet<string>();
            var agent = new Algorithms.Random(settings);
            for (var t = 0; t < NUM_TRAIL; t++)
            {
                Console.WriteLine("Trail: {0}", t);
                var cand = Seed.copy();

                for (var step = 0; step < TOTAL_RULE ; step++)
                {
                    var opt = agent.ChooseOption(cand);
                    if (opt == null)
                    {
                        Console.WriteLine("Fail on step {0}", step+1);
                        return;
                    }

                    agent.ApplyOption(opt, cand, true);
                }
                var carboxOpt = agent.ChooseCarboxOption(cand);
                if (carboxOpt == null)
                {
                    Console.WriteLine("Fail on finding final carbox");
                    return;
                }
                agent.ApplyOption(carboxOpt, cand, true);
                var candSmile = OBFunctions.moltoSMILES(OBFunctions.designgraphtomol(cand.graph));
                var linkerName = AbstractAlgorithm.GetLinkerName(cand);
                Console.WriteLine(candSmile);
                Console.WriteLine(linkerName);
                if (linkerSet.Contains(linkerName))
                {
                    t--;
                    continue;
                }
                linkerSet.Add(linkerName);
                var coeff = Path.Combine(_runDirectory, "data", "linker" + linkerName + ".coeff");
                var lmpdat = Path.Combine(_runDirectory, "data", "linker" + linkerName + ".lmpdat");
                agent.Converter.moltoUFF(OBFunctions.designgraphtomol(cand.graph), coeff, lmpdat, false, 100);
                server.CalculateProperty(linkerName);

                //mutex.WaitOne();
                jobBuffer.Add(linkerName, AbstractAlgorithm.Rand.NextDouble());
                //mutex.ReleaseMutex();
            }
            jobBuffer.Add("finish", 1.0);
            
            
        }
        
        public override string text => "RandomTrail Search Runner";
        
    }
}