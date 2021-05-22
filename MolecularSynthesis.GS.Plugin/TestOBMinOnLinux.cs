using System.Collections.Generic;
using System.Threading;
using GraphSynth;
using GraphSynth.Representation;
using GraphSynth.Search;
using Priority_Queue;
using System;
using System.Collections;
using System.Security.Cryptography.X509Certificates;
using System.IO;
using MolecularSynthesis.GS.Plugin;
using System.Linq;
using OpenBabel;
using OpenBabelFunctions;
using MolecularSynthesis.GS.Plugin;
using System.Diagnostics;

namespace MolecularSynthesis.GS.Plugin
{
    public class TestOBMinOnLinux : SearchProcess
    {
        // give desiredMoment
        // [] desiredMoment = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        static double[] desiredLenghtAndRadius = new double[] { 300, 50 };
        static Random rnd = new Random(0);

        public TestOBMinOnLinux(GlobalSettings settings) : base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 2;
            AutoPlay = true;
        }


        public override string text
        {
            get { return "TestOBMinOnLinux"; }
        }

        protected override void Run()
        {

            //TreeCandidate StartState = new TreeCandidate(seedCandidate);

            //StartState.S = 0;
            //StartState.n = 0;
            //StartState.UCB = double.MaxValue;
            //StartState.Children = new List<TreeCandidate>();

            //var option0 = rulesets[0].recognize(StartState.graph);
            //var option1 = rulesets[1].recognize(StartState.graph);
            //var option2 = rulesets[2].recognize(StartState.graph);

            ////option0[6].apply(StartState.graph, null);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[6].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[6]);

            // add a stopwatch to record time
            var timer = new Stopwatch();
            timer.Start();

            TreeCandidate StartState = new TreeCandidate(seedCandidate);
            var candidate = (TreeCandidate)StartState.copy();

            var option0 = rulesets[0].recognize(candidate.graph);
            var option1 = rulesets[1].recognize(candidate.graph);
            var option2 = rulesets[2].recognize(candidate.graph);

            option0 = rulesets[0].recognize(candidate.graph);
            option0[5].apply(candidate.graph, null);

            //option0 = rulesets[0].recognize(candidate.graph);
            //option0[2].apply(candidate.graph, null);

            //option0 = rulesets[0].recognize(candidate.graph);
            //option0[6].apply(candidate.graph, null);

            //option1 = rulesets[1].recognize(candidate.graph);
            //option1[1].apply(candidate.graph, null);

            option2 = rulesets[2].recognize(candidate.graph);
            option2[0].apply(candidate.graph, null);

            double PoreSize = Evaluation.GetPoreSize(candidate.graph);
            Console.WriteLine("PoreSize= "+ PoreSize);

            //var resultMol = OBFunctions.designgraphtomol(candidate.graph);
            //var conv = new OBConversion();
            //conv.SetInAndOutFormats("pdb", "mol");
            //Console.WriteLine("writing TestOBD.mol");
            //conv.WriteFile(resultMol, Path.Combine("/nfs/hpc/share/zhangho2/MolecularSynthesis/output", "TestOBD.mol"));







            timer.Stop();
            TimeSpan ts = timer.Elapsed;
            //string foo = "Time taken: " + timeTaken.ToString(@"m\:ss\.fff");

            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
            ts.Hours, ts.Minutes, ts.Seconds,
            ts.Milliseconds / 10);
            Console.WriteLine("RunTime " + elapsedTime);


        }

    }

}






