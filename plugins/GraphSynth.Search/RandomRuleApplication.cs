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
        private Computation computation;
        private StreamWriter writer;
        private LearningServer server;
        private MessageClient client;


        private const int NUM_EPOCH = 2;
        private const int NUM_TRAIL = 2;
        private const int TOTAL_RULE_MIN = 6;
        private const int TOTAL_RULE_MAX = 16;
        private const string CARBOXTYPE = "estimator";
        
        //private static Mutex mutex = new Mutex();


        public RandomRuleApplication(GlobalSettings settings): base(settings) 
        {
            RequireSeed = true;
            RequiredNumRuleSets = 1;
            AutoPlay = true;
            
            Seed = new candidate(OBFunctions.tagconvexhullpoints(settings.seed), settings.numOfRuleSets);
            
            _runDirectory = Path.Combine(settings.OutputDirAbs, "RandomRuleApplication", CARBOXTYPE);
            
            if (Directory.Exists(_runDirectory))
                Directory.Delete(_runDirectory, true);
            Directory.CreateDirectory(_runDirectory);

            jobBuffer = new JobBuffer(_runDirectory);
            
            var learnDirectory = Path.Combine(settings.OutputDirAbs, "morfLearn");
            computation = new Computation(_runDirectory, learnDirectory, "point", "stiff");
            writer = new StreamWriter(Path.Combine(_runDirectory, CARBOXTYPE + ".txt"));

            System.Random rnd = new System.Random();

            var port = rnd.Next(1, 65535);
            port = 9999;
            server = new LearningServer(learnDirectory, port, _runDirectory);
            client = new MessageClient(port);
        }

        protected override void Run()
        {
            //server.StartOnlineServer();

            client.Connect();
            //client.SendMessage("[Time]");
            //Thread.Sleep(5000);


            Thread generateLinkers = new Thread(Generate);
            Thread autoReleaseBuffer = new Thread(AutoSubmitSimulation);
            
            generateLinkers.Start();
            autoReleaseBuffer.Start();

            generateLinkers.Join();
            autoReleaseBuffer.Join();

            client.DisConnect();
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
                var allFinished = jobBuffer.Check_finised(computation, writer, client);
                if (allFinished && allSubmitted)
                {
                    writer.Close();
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
                        candidate cand = null;
                        while(cand == null)
                        {
                            cand = Seed.copy();
                            Console.WriteLine("Total Intermediate Rules: {0}", total_rule);
                            for (var step = 0; step < total_rule; step++)
                            {
                                cand = agent.ChooseAndApplyOption(cand);
                                if (cand == null)
                                {
                                    Console.WriteLine("Fail on step {0}", step+1);
                                    break;
                                }
                            }
                            if (cand == null)
                                continue;
                            //var carboxOpt = agent.ChooseAndApplyOption(cand);
                            //var carboxOpt = agent.ChooseAndApplyCarboxOptionBestAngle(cand);
                            cand = agent.ChooseAndApplyCarboxOptionUsingEstimator(cand, computation, client, _runDirectory);
                            if (cand == null)
                                Console.WriteLine("Fail on finding final carbox");
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

                        double piority = 0;
                        if (CARBOXTYPE == "estimator")
                        {
                            piority = - Convert.ToDouble(client.SendMessage("[Predict]" + " " + linkerName));
                        }
                        else
                        {
                            piority = AbstractAlgorithm.Rand.NextDouble();
                            computation.CalculateFeature(linkerName);
                        }

                        //mutex.WaitOne();
                        jobBuffer.Add(linkerName, piority, e);
                        //mutex.ReleaseMutex();
                    }
                }
            }
            jobBuffer.AllSubmitFlag = true;
        }

        private void GenerateFixed()
        {
            var agent = new Algorithms.Random(settings);
            var linkerBeforeCarboxSet = new HashSet<string>();
            for (var t = 0; t < NUM_TRAIL; t++)
            {
                Console.WriteLine("Trail: {0}", t);
                for (var total_rule = TOTAL_RULE_MIN; total_rule < TOTAL_RULE_MAX; total_rule++)
                {
                    candidate cand = null;
                    while(cand == null)
                    {
                        cand = Seed.copy();
                        Console.WriteLine("Total Intermediate Rules: {0}", total_rule);
                        for (var step = 0; step < total_rule; step++)
                        {
                            cand = agent.ChooseAndApplyOption(cand);
                            if (cand == null)
                            {
                                Console.WriteLine("Fail on step {0}", step+1);
                                break;
                            }
                        }
                    }
                    var linkerName = AbstractAlgorithm.GetLinkerName(cand);
                    Console.WriteLine(linkerName);
                    if (linkerBeforeCarboxSet.Contains(linkerName))
                    {
                        total_rule--;
                        continue;
                    }
                    linkerBeforeCarboxSet.Add(linkerName);
                }
            }
            Console.WriteLine(linkerBeforeCarboxSet.Count);
        }
        
        public override string text => "RandomTrail Search Runner";
        
    }
}