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

        private bool allSubmitFlag;


        private const int NUM_EPOCH = 10;
        private const int NUM_TRAIL = 10;
        private const int TOTAL_RULE_MIN = 6;
        private const int TOTAL_RULE_MAX = 15;
        private const string CARBOXTYPE = "estimator";
        
        //private static Mutex sendMessageMutex = new Mutex();


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
            //port = 9999;
            server = new LearningServer(learnDirectory, port, _runDirectory);
            client = new MessageClient(port);
        }

        protected override void Run()
        {
            server.StartOnlineServer();
            client.Connect(10);

            //Thread generateLinkers = new Thread(Generate);
            Thread generateLinkers = new Thread(GenerateFixed);
            Thread autoReleaseBuffer = new Thread(AutoSubmitSimulation);
            
            generateLinkers.Start();
            autoReleaseBuffer.Start();

            generateLinkers.Join();
            autoReleaseBuffer.Join();

            client.DisConnect();
            server.ShutDownOnlineServer();
        }
        
        private void AutoSubmitSimulation()
        {
            while (true)
            {
                //mutex.WaitOne();
                if (jobBuffer.CanFeedIn())
                {
                    jobBuffer.Simulate();
                }
                var allFinished = jobBuffer.Check_finised(computation, writer, client);
                if (allFinished && allSubmitFlag)
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
                    for (var total_rule = TOTAL_RULE_MIN; total_rule < TOTAL_RULE_MAX + 1; total_rule++)
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
                if (CARBOXTYPE == "estimator")
                {
                    while(true)
                    {
                        var on_simulation = jobBuffer.Num_simulating();
                        if (on_simulation == 0)
                            break;
                        Console.WriteLine("Wait for current {0} linkers to finish simulation....", on_simulation);
                        Thread.Sleep(10000);
                    }
                    client.SendMessage("[FitModel]");
                }
            }
            allSubmitFlag = true;
        }

        private void GenerateFixed()
        {
            var agent = new Algorithms.Random(settings);
            var linkerBeforeCarboxDict = new Dictionary<string, candidate>();
            for (var t = 0; t < NUM_TRAIL; t++)
            {
                Console.WriteLine("Trail: {0}", t);
                for (var total_rule = TOTAL_RULE_MIN; total_rule < TOTAL_RULE_MAX + 1; total_rule++)
                {
                    candidate cand = null;
                    candidate finalCand = null;
                    while(cand == null || finalCand == null)
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
                        finalCand = agent.ChooseAndApplyCarboxOption(cand);
                        if (finalCand == null)
                            Console.WriteLine("Fail on finding final carbox");
                    }
                    var linkerName = AbstractAlgorithm.GetLinkerName(cand);
                    Console.WriteLine(linkerName);
                    if (linkerBeforeCarboxDict.ContainsKey(linkerName))
                    {
                        total_rule--;
                        continue;
                    }
                    linkerBeforeCarboxDict[linkerName] = cand;
                }
            }
            Console.WriteLine(linkerBeforeCarboxDict.Count);
            for (var e = 0; e < NUM_EPOCH; e++)
            {
                Console.WriteLine("Epoch: {0}", e);
                foreach(var item in linkerBeforeCarboxDict)
                {
                    var submitCand = agent.ChooseAndApplyCarboxOptionUsingEstimator(item.Value, computation, client, _runDirectory, e);
                    //cand = agent.ChooseAndApplyCarboxOptionBestAngle(item.Value());
                    //cand = agent.ChooseAndApplyCarboxOptionUsingEstimator(item.Value(), computation, client, _runDirectory);
                    if (submitCand == null)
                    {
                        Console.WriteLine("Fail on finding final carbox, should never happen");
                        Environment.Exit(0);
                    }
                    var linkerName = AbstractAlgorithm.GetLinkerName(submitCand) + "-E" + e.ToString();
                    var coeff = Path.Combine(_runDirectory, "data", "linker" + linkerName + ".coeff");
                    var lmpdat = Path.Combine(_runDirectory, "data", "linker" + linkerName + ".lmpdat");
                    agent.Converter.moltoUFF(OBFunctions.designgraphtomol(submitCand.graph), coeff, lmpdat, false, 100);

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

            allSubmitFlag = true;
        }
        
        public override string text => "RandomTrail Search Runner";
        
    }
}