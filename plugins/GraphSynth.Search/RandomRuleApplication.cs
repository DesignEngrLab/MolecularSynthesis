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
        private StreamWriter sw;

        private const int NUM_EPOCH = 10;
        private const int NUM_TRAIL = 10;
        private const int TOTAL_RULE_MIN = 6;
        private const int TOTAL_RULE_MAX = 16;
        
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
            
            if (Directory.Exists(_runDirectory))
                Directory.Delete(_runDirectory, true);
            Directory.CreateDirectory(_runDirectory);

            jobBuffer = new JobBuffer(_runDirectory);
            
            var learnDirectory = Path.Combine(settings.OutputDirAbs, "morfLearn");
            server = new LearningServer(_runDirectory, learnDirectory, "point", "stiff");
            sw = new StreamWriter(Path.Combine(_runDirectory, "stiff" + ".txt"));
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
                var allFinished = jobBuffer.Check_finised(server, sw);
                if (allFinished && allSubmitted)
                {
                    sw.Close();
                    break;
                }

                //mutex.ReleaseMutex();
            }
        }

        private void Generate()
        {
            var agent = new Algorithms.Random(settings);
            var linkerSet = new HashSet<string>();
            for (var e = 0; e < NUM_EPOCH; e++)
            {
                Console.WriteLine("Epoch: {0}", e);
                for (var t = 0; t < NUM_TRAIL; t++)
                {
                    Console.WriteLine("Trail: {0}", t);
                    for (var total_rule = TOTAL_RULE_MIN; total_rule < TOTAL_RULE_MAX; total_rule++)
                    {
                        var cand = Seed.copy();
                        while(true)
                        {
                            var successFlag = true;
                            Console.WriteLine("Total Intermediate Rules: {0}", total_rule);
                            for (var step = 0; step < total_rule; step++)
                            {
                                var opt = agent.ChooseOption(cand);
                                if (opt == null)
                                {
                                    Console.WriteLine("Fail on step {0}", step+1);
                                    successFlag = false;
                                    cand = Seed.copy();
                                    break;
                                }
                                agent.ApplyOption(opt, cand, true);
                            }
                            if (successFlag == false)
                                continue;
                            var carboxOpt = agent.ChooseCarboxOption(cand);
                            if (carboxOpt == null)
                            {
                                Console.WriteLine("Fail on finding final carbox");
                                successFlag = false;
                                cand = Seed.copy();
                            }
                            if (successFlag == false)
                                continue;
                            agent.ApplyOption(carboxOpt, cand, true);
                            break;
                        }
                        
                        var candSmile = OBFunctions.moltoSMILES(OBFunctions.designgraphtomol(cand.graph));
                        var linkerName = AbstractAlgorithm.GetLinkerName(cand);
                        //Console.WriteLine(candSmile);
                        Console.WriteLine(linkerName);
                        if (linkerSet.Contains(linkerName))
                        {
                            total_rule--;
                            continue;
                        }
                        linkerSet.Add(linkerName);
                        var coeff = Path.Combine(_runDirectory, "data", "linker" + linkerName + ".coeff");
                        var lmpdat = Path.Combine(_runDirectory, "data", "linker" + linkerName + ".lmpdat");
                        agent.Converter.moltoUFF(OBFunctions.designgraphtomol(cand.graph), coeff, lmpdat, false, 100);
                        server.CalculateFeature(linkerName);

                        //mutex.WaitOne();
                        jobBuffer.Add(linkerName, AbstractAlgorithm.Rand.NextDouble(), e);
                        //mutex.ReleaseMutex();
                    }
                }
            }
            jobBuffer.AllFinishFlag = true;
        }
        
        public override string text => "RandomTrail Search Runner";
        
    }
}