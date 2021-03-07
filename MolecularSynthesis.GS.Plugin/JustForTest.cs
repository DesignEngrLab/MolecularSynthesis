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

using System.Threading;
using System.Threading.Tasks;

namespace MolecularSynthesis.GS.Plugin
{
    public class JustForTest : SearchProcess
    {
        // give desiredMoment
        // [] desiredMoment = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        static double[] desiredLenghtAndRadius = new double[] { 300, 50 };
        //static Random rnd = new Random(0);

        public JustForTest(GlobalSettings settings) : base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 2;
            AutoPlay = true;
        }


        public override string text
        {
            get { return "JustForTest"; }
        }

        protected override void Run()
        {

            TreeCandidate StartState = new TreeCandidate(seedCandidate);

            StartState.S = 0;
            StartState.n = 0;
            StartState.UCB = double.MaxValue;
            StartState.Children = new List<TreeCandidate>();            

            // add a stopwatch to record time
            var timer = new Stopwatch();
            timer.Start();

            // Randomly generate .mol and .xyz files
            int TotalNumber = 8;
            var rand = new Random();

            //TreeCandidate StartState = new TreeCandidate(seedCandidate);
            // startStateCopies = new TreeCandidate[TotalNumber];
            //for (int i = 0; i < TotalNumber; i++)
            //    startStateCopies[i] = (TreeCandidate)StartState.copy();
            
            Parallel.For(0, TotalNumber, i =>
            //for (int i = 0; i < TotalNumber; i++)
            {
                //var candidate = startStateCopies[i];
                var candidate = (TreeCandidate)StartState.copy();

                var option0 = rulesets[0].recognize(candidate.graph,false);
                var option1 = rulesets[1].recognize(candidate.graph, false);
                var option2 = rulesets[2].recognize(candidate.graph, false);

                option0 = rulesets[0].recognize(candidate.graph,false);
                option0[6].apply(candidate.graph, null);
                //candidate.addToRecipe(option0[6]);
                Debug.WriteLine("\r\nInitialnizing"+ Thread.CurrentThread.ManagedThreadId);

                //foreach (var rule in candidate.recipe)
                //{
                //    Debug.WriteLine("ThreadNumber" + Thread.CurrentThread.ManagedThreadId + " "+rule.ruleSetIndex+" "+ rule.ruleNumber);
                //    Debug.WriteLine("------BezeneRing");
                    
                //}
                


                for (int j = 0; j < 4; j++)
                {
                    //rnd.Next(0, 2); 0 or 1
                    var RuleSetNumber = rand.Next(0, 2);
                    var TotalOption = rulesets[RuleSetNumber].recognize(candidate.graph, false).Count;
                    var OptionNumber = rand.Next(0, TotalOption);
                    rulesets[RuleSetNumber].recognize(candidate.graph,false)[OptionNumber].apply(candidate.graph, null);
                    Debug.WriteLine("growing" + Thread.CurrentThread.ManagedThreadId);
                    //candidate.addToRecipe(rulesets[RuleSetNumber].recognize(candidate.graph,false)[OptionNumber]);

                    //foreach (var rule in candidate.recipe)
                    //{
                    //    Debug.WriteLine("ThreadNumber" + Thread.CurrentThread.ManagedThreadId + " " + rule.ruleSetIndex + " " + rule.ruleNumber);
                    //    Debug.WriteLine("------OtherRules");
                    //}


                }
                
                option2 = rulesets[2].recognize(candidate.graph, false);
                option2[0].apply(candidate.graph, null);
                //candidate.addToRecipe(option2[0]);
                Debug.WriteLine("Stop growing" + Thread.CurrentThread.ManagedThreadId);


                //foreach (var rule in candidate.recipe)
                //{
                //    Debug.WriteLine("ThreadNumber" + Thread.CurrentThread.ManagedThreadId + " " + rule.ruleSetIndex + " " + rule.ruleNumber);
                //    Debug.WriteLine("------Carboxylate");
                //}

                var resultMol = OBFunctions.designgraphtomol(candidate.graph);
                Debug.WriteLine("complete graph to mol" + Thread.CurrentThread.ManagedThreadId);
                resultMol = OBFunctions.InterStepMinimize(resultMol);
                Debug.WriteLine("start updating" + Thread.CurrentThread.ManagedThreadId);
                OBFunctions.updatepositions(candidate.graph, resultMol);

                var FinalResultMol = OBFunctions.designgraphtomol(candidate.graph);

                var conv = new OBConversion();
                conv.SetInAndOutFormats("pdb", "mol");

                string name = ".mol";

                Debug.WriteLine("start writing .mol file" + Thread.CurrentThread.ManagedThreadId);
                name = Convert.ToString(Thread.CurrentThread.ManagedThreadId) + name;
                conv.WriteFile(FinalResultMol, Path.Combine("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\examples", name));
                       



                //string name2 = ".xyz";
                //name2 = Convert.ToString(i) + name2;

                //using (Process proc = new Process())
                //{
                //    //"C:\Program Files\OpenBabel-3.1.1\obabel.exe"
                //    proc.StartInfo.FileName = "C:\\Program Files\\OpenBabel-3.1.1\\obabel.exe";
                //    proc.StartInfo.Arguments = name + " -O " + name2;
                //    proc.StartInfo.WorkingDirectory = "C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\examples";

                //    //proc.StartInfo.RedirectStandardError = true;
                //    //proc.StartInfo.UseShellExecute = false;
                //    proc.StartInfo.RedirectStandardOutput = true;
                //    //proc.StartInfo.RedirectStandardInput = false;

                //    Console.Write("starting Convert...");
                //    proc.Start();

                //    //minimizeOutput = proc.StandardOutput.ReadToEnd();
                //    proc.WaitForExit();
                //}
            });
            //}
            timer.Stop();
            TimeSpan ts = timer.Elapsed;
            //string foo = "Time taken: " + timeTaken.ToString(@"m\:ss\.fff");

            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
            ts.Hours, ts.Minutes, ts.Seconds,
            ts.Milliseconds / 10);
            Console.WriteLine("RunTime " + elapsedTime);

            //SearchIO.output("Total running time:" , elapsedTime);
            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[0].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[0]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[1].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[1]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[2].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[2]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[3].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[3]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[4].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[4]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[5].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[5]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[6].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[6]);

            //// 3(1),7(2),10(3),12(4),16(5),22(6),*26(7),*31(8),*35(9)
            //option1 = rulesets[1].recognize(StartState.graph);
            //option1[3].apply(StartState.graph, null);
            //StartState.addToRecipe(option1[3]);

            //option1 = rulesets[1].recognize(StartState.graph);
            //option1[25].apply(StartState.graph, null);
            //StartState.addToRecipe(option1[25]);

            //option2 = rulesets[2].recognize(StartState.graph);
            //option2[0].apply(StartState.graph, null);
            //StartState.addToRecipe(option2[0]);

            ////Save("XYZ.gxml",StartState.graph);



            //var resultMol = OBFunctions.designgraphtomol(StartState.graph);
            //resultMol = OBFunctions.InterStepMinimize(resultMol);
            //OBFunctions.updatepositions(StartState.graph, resultMol);
            //var FinalResultMol= OBFunctions.designgraphtomol(StartState.graph);
            //var conv = new OBConversion();
            //conv.SetInAndOutFormats("pdb", "xyz");
            //conv.WriteFile(FinalResultMol, Path.Combine("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\output", "Test112.xyz"));

            //var score = Evaluation.distance(StartState, desiredLenghtAndRadius);

            //double TotalMass = Evaluation.TotalAtomMass(StartState);

            //SearchIO.output(score + "  HAHA");

            //for (int j = 0; j < StartState.recipe.Count; j++)
            //{
            //    SearchIO.output(StartState.recipe[j].ruleSetIndex + " " + StartState.recipe[j].optionNumber);

            //}

            // test for Rule from 1-6 from RS0 , no problem with this part
            //for (int i = 0; i < 6; i++)
            //{
            //    TreeCandidate candidate = StartState;

            //    option0[i].apply(candidate.graph, null);
            //    //option1[].apply(candidate.graph, null);
            //    option2[0].apply(candidate.graph, null);

            //    var resultMol0 = OBFunctions.designgraphtomol(candidate.graph);
            //    resultMol0 = OBFunctions.InterStepMinimize(resultMol0);

            //    OBFunctions.updatepositions(seedGraph, resultMol0);

            //    var score0 = Evaluation.distance(candidate, desiredLenghtAndRadius);
            //    SearchIO.output(score0+"RS0");
            //}


            // test for Rule from 1-9 from RS1 , no problem with this part except rule 9
            //option0[6].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[6]);

            //option2 = rulesets[2].recognize(StartState.graph);
            //var option1 = rulesets[1].recognize(StartState.graph);

            //for (int i = 0; i < 8; i++)
            //{
            //    TreeCandidate candidate = StartState;

            //    option1[i*8].apply(candidate.graph, null);
            //    option2[0].apply(candidate.graph, null);

            //    var resultMol1 = OBFunctions.designgraphtomol(candidate.graph);
            //    resultMol1 = OBFunctions.InterStepMinimize(resultMol1);

            //    OBFunctions.updatepositions(seedGraph, resultMol1);

            //    var score1 = Evaluation.distance(candidate, desiredLenghtAndRadius);
            //    SearchIO.output(score1 + "RS1");

            //}

            //option1[32].apply(StartState.graph,null);
            //option1[8].apply(StartState.graph, null);
            //option1[16].apply(StartState.graph, null);
            //option1[64].apply(StartState.graph, null);
            //option2[0].apply(StartState.graph, null);
            //var resultMol = OBFunctions.designgraphtomol(StartState.graph);
            //resultMol = OBFunctions.InterStepMinimize(resultMol);

            //OBFunctions.updatepositions(seedGraph, resultMol);

            //var score = Evaluation.distance(StartState, desiredLenghtAndRadius);
            //SearchIO.output(score + "HAHA");

            //for (int i = 0; i < iteration; i++)
            //{


            //TreeCandidate current = StartState;
            //current.S = Rollout(current);
            //SearchIO.output("distance:" + current.S);

            //}
        }






    }

}






