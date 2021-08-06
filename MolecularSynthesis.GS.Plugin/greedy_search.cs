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
using System.Threading.Tasks;
using System.Threading;

namespace MolecularSynthesis.GS.Plugin
{
    public class greedy_search : SearchProcess
    {
        // give desiredMoment
        // [] desiredMoment = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        //static double[] desiredLenghtAndRadius = new double[] { 300, 50 };
        static Random rnd = new Random(0);
        //static double[] desiredLenghtAndRadius = new double[] { 245.277, 89.53 };
        // RS0 R7 R3; RS1 R1 R2

        static double[] desiredLenghtAndRadius = new double[] { 658.15, 94.29 };
        // RS0 R3 R4 R5 R6

        public greedy_search(GlobalSettings settings) : base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 2;
            AutoPlay = true;
        }

        static TreeCandidate noneparallel = new TreeCandidate(new candidate());

        public override string text
        {
            get { return "greedy_search"; }
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

            // Randomly generate .mol and .xyz files


            // thread number
            int TotalNumber = 12;
            var rand = new Random();

            //TreeCandidate StartState = new TreeCandidate(seedCandidate);

            TreeCandidate StartState = new TreeCandidate(seedCandidate);

            StartState.S = 0;
            StartState.n = 0;
            StartState.UCB = double.MaxValue;
            StartState.Children = new List<TreeCandidate>();
            StartState.Parent = null;

            List<string> resultCollector = new List<string>();


            // thread number on HPC cn-h-5(5,6,7,8) is 58, when apply for 44 core



            // greedy search
            // 1.make an empty list to operate queue option
            // 2.start from seed, generate children, then do evaluation and put them into the empty list
            // 3.sort list, get the biggest one, then expand this candidate
            // 4.do until reach maximim iteration or find the target

            // parallel greedy search(acutually beam search)
            // 1.make an empty list to operate queue option
            // 2.start from seed, generate children, then do evaluation and put them into the empty list
            // 3.sort list, each thread will get the biggest one serially, then expand those candidates and put their children into the list parallel
            // 4.do until reach maximim iteration or find the target

            // Epsilon-Greedy search

            //for (int i = 0; i < TotalNumber; i++)
            //{

            List<TreeCandidate> greedySearchQueue = new List<TreeCandidate>();

            //try
            //{



            //Parallel.For(0, TotalNumber, count =>
            //{

            var candidate = (TreeCandidate)StartState.copy();
            //var bestSolution = candidate;

            //greedySearchQueue.Add(candidate);

            var option0 = rulesets[0].recognize(candidate.graph);
            //var option1 = rulesets[1].recognize(candidate.graph);
            //var option2 = rulesets[2].recognize(candidate.graph);

            // create the initial list for the next expansion

            for (int i = 0; i < option0.Count; i++)
            {
                

                option0= rulesets[0].recognize(candidate.graph);
                option0[i].apply(candidate.graph, null);
                candidate.addToRecipe(option0[i]);

                candidate.Parent = StartState;

                var ForEvaluationOnly= (TreeCandidate)candidate.copy();
                candidate.S= simpleEvaluation(ForEvaluationOnly);


                greedySearchQueue.Add(candidate);

                candidate = (TreeCandidate)StartState.copy();
            }




            //foreach (var option in option0)
            //{
            //    var child = (TreeCandidate)candidate.copy();

            //    option0 = rulesets[0].recognize(child.graph);
            //    option2 = rulesets[2].recognize(child.graph);

            //    option.apply(child.graph, null);
            //    child.addToRecipe(option);

            //    option2 = rulesets[2].recognize(child.graph);
            //    option2[0].apply(child.graph, null);
            //    candidate.addToRecipe(option0[2]);

            //    var resultMol = OBFunctions.designgraphtomol(candidate.graph);
            //    resultMol = OBFunctions.InterStepMinimize(resultMol);
            //    OBFunctions.updatepositions(candidate.graph, resultMol);




            //    child.Parent = candidate;

            //    child.S = simpleEvaluation(child);

            //    candidate.Children.Add(child);
            //    greedySearchQueue.Add(child);

            //}

            Console.WriteLine("-----------------before sort-----------------");
            foreach (var child in greedySearchQueue)
            {
                Console.WriteLine(child.S);
            }

            // sort list based on S value
            SortByS sortByName = new SortByS();
            greedySearchQueue.Sort(sortByName);

            Console.WriteLine("--------------after sort----------------------");
            foreach (var child in greedySearchQueue)
            {
                Console.WriteLine(child.S);
            }

            Console.WriteLine("-----------------------------------------");
            //var QueueInOrder = greedySearchQueue.OrderBy(TreeCandidate => TreeCandidate.S);
            //greedySearchQueue.Sort();
            //greedySearchQueue.Remove(candidate);

            List<TreeCandidate> waitingQueue = new List<TreeCandidate>();

            //Console.WriteLine(GetSuccessor(greedySearchQueue[0]).Count);



            //for (int i = 0; i < greedySearchQueue.Count; i++)
            //{

            //    var conv = new OBConversion();
            //    lock (noneparallel)
            //        conv.SetInAndOutFormats("pdb", "mol");

            //    //int ThreadNumber = System.Threading.Thread.CurrentThread.ManagedThreadId;
            //    string filename = "HHH" + i.ToString() + ".mol";

            //    var HHH = OBFunctions.designgraphtomol(greedySearchQueue[i].graph);
            //    lock (noneparallel)
            //        conv.WriteFile(HHH, Path.Combine("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\examples", filename));

            //}


            // greedy search begin
            for (int i = 0; i < 2; i++)
            {
                // determine the total number of candidates for 40 threads to evaluate
                var optionCounter = GetSuccessor(greedySearchQueue[0]).Count;
                var candidateCounter = 1;
                //waitingQueue.Add(greedySearchQueue[0]);


                while (optionCounter <= TotalNumber)
                {
                    var nextCandidateOptionNumber = GetSuccessor(greedySearchQueue[candidateCounter]).Count;

                    optionCounter = optionCounter + nextCandidateOptionNumber;
                    candidateCounter = candidateCounter + 1;

                    if (optionCounter >= TotalNumber)
                    {
                        optionCounter = optionCounter - nextCandidateOptionNumber;
                        candidateCounter = candidateCounter - 1;
                        break;
                    }

                }

                // build the list for parallel evaluation, first candidate ,second candidate, third candidate ...
                // and update queue
                for (int j = 0; j < candidateCounter; j++)
                {


                    for (int k = 0; k < GetSuccessor(greedySearchQueue[j]).Count; k++)
                    {
                        //greedySearchQueue.Add(GetSuccessor(greedySearchQueue[j])[k]);
                        waitingQueue.Add(GetSuccessor(greedySearchQueue[j])[k]);
                        greedySearchQueue[j].Children.Add(GetSuccessor(greedySearchQueue[j])[k]);

                    }

                }

                // start parallel evaluation
                Parallel.For(0, waitingQueue.Count,i=> { waitingQueue[i].S = simpleEvaluation(waitingQueue[i]); });

                // put new candidates into greedy search queue and clean the waiting queue 

                for (int j = 0; j < waitingQueue.Count; j++)
                {
                    greedySearchQueue.Add(waitingQueue[j]);
                }

                waitingQueue.Clear();

                // remove the older candidate
                for (int m = 0; m < candidateCounter; m++)
                {
                    greedySearchQueue.RemoveAt(m);
                }

                // sort list for next iteration
                greedySearchQueue.Sort(sortByName);



            }

            foreach (var child in greedySearchQueue)
            {
                Console.WriteLine(child.S);
            }





            //for (int i = 0; i < 10; i++)
            //{



            //    //greedySearchQueue.Add();
            //    //greedySearchQueue.Sort();



            //    // generate backbone rules
            //    var BBnumber = rand.Next(1, 11);

            //    for (int j = 0; j < BBnumber; j++)
            //    {
            //        option0 = rulesets[0].recognize(candidate.graph);
            //        var totalOptionNumber = option0.Count();
            //        var OptionNumber = rand.Next(0, totalOptionNumber);
            //        option0[OptionNumber].apply(candidate.graph, null);
            //        candidate.addToRecipe(option0[OptionNumber]);
            //    }


            //    // generate functional group rules

            //    //option1 = rulesets[1].recognize(candidate.graph);
            //    //var PossibleFGnumber = option1.Count;
            //    //var totalFGnumber = rand.Next(0, PossibleFGnumber);


            //    //for (int j = 0; j < totalFGnumber; j++)
            //    //{
            //    //    var OptionNumber2= rand.Next(0, PossibleFGnumber);
            //    //    option1[OptionNumber2].apply(candidate.graph,null);
            //    //    candidate.addToRecipe(option1[OptionNumber2]);


            //    //}                                

            //    option2 = rulesets[2].recognize(candidate.graph);
            //    option2[0].apply(candidate.graph, null);
            //    //candidate.addToRecipe(option2[0]);

            //    var resultMol = OBFunctions.designgraphtomol(candidate.graph);
            //    resultMol = justMinimize(resultMol);
            //    OBFunctions.updatepositions(candidate.graph, resultMol);

            //    var score = Evaluation.distance(candidate, desiredLenghtAndRadius);
            //    resultCollector.Add(score.ToString());

            //    //----------------------------------------
            //    //var score = Evaluation.TotalAtomMass(candidate);
            //    //resultCollector.Add(score.ToString());


            //    // check the candidate in .mol

            //    //var FinalResultMol = OBFunctions.designgraphtomol(candidate.graph);

            //    //var conv = new OBConversion();
            //    //lock (noneparallel)
            //    //    conv.SetInAndOutFormats("pdb", "mol");

            //    //string name = ".mol";
            //    //name = Convert.ToString(System.Threading.Thread.CurrentThread.ManagedThreadId) + name;
            //    //lock (noneparallel)
            //    //    conv.WriteFile(FinalResultMol, Path.Combine("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\examples", name));

            //}
            //});

            //}


            //catch (Exception exc)
            //{
            //    Console.WriteLine("exc.message");
            //}

            timer.Stop();
            TimeSpan ts = timer.Elapsed;
            //string foo = "Time taken: " + timeTaken.ToString(@"m\:ss\.fff");

            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
            ts.Hours, ts.Minutes, ts.Seconds,
            ts.Milliseconds / 10);
            Console.WriteLine("RunTime " + elapsedTime);

            System.IO.File.WriteAllLines(@"C:\Users\zhang\source\repos\MolecularSynthesis\examples\RandomSearchRecord.txt", resultCollector);


            //}
        }

        public class SortByS : IComparer<TreeCandidate>
        {
            public int Compare(TreeCandidate x, TreeCandidate y)
            {
                return x.S.CompareTo(y.S);
            }
        }


        public List<TreeCandidate> GetSuccessor(TreeCandidate current)
        {
               
            

            List<TreeCandidate> candidateList = new List<TreeCandidate>();
            var option0 = rulesets[0].recognize(current.graph);
            int PotenialOptionNumber = option0.Count;


            for (int i = 0; i < PotenialOptionNumber; i++)
            {
                var child = (TreeCandidate)current.copy();
                child.Parent = current;
                child.Children = new List<TreeCandidate>();
                child.n = 0;
                child.S = 0;
                child.UCB = double.MaxValue;

                option0 = rulesets[0].recognize(child.graph);
                option0[i].apply(child.graph, null);
                child.addToRecipe(option0[i]);

                candidateList.Add(child);
                current.Children.Add(child);
            }

                     

           


            return candidateList;
        }


        public double simpleEvaluation(TreeCandidate candidate)
        {
            var option0 = rulesets[0].recognize(candidate.graph);
            var option1 = rulesets[1].recognize(candidate.graph);
            var option2 = rulesets[2].recognize(candidate.graph);

            //var option2 = rulesets[2].recognize(candidate.graph);

            option2[0].apply(candidate.graph, null);
            candidate.addToRecipe(option2[0]);

            var resultMol = OBFunctions.designgraphtomol(candidate.graph);
            resultMol = justMinimize(resultMol);
            OBFunctions.updatepositions(candidate.graph, resultMol);

            var score = Evaluation.distance(candidate, desiredLenghtAndRadius);
            //resultCollector.Add(score.ToString());

            //candidate.undoLastRule();

            return 10000 - score;
        }



        private static OBMol justMinimize(OBMol mol)
        {
            var conv = new OBConversion();
            lock (noneparallel)
                conv.SetInAndOutFormats("pdb", "mol");

            //lock (noneparallel) 

            int ThreadNumber = System.Threading.Thread.CurrentThread.ManagedThreadId;
            string filename = "Test" + ThreadNumber.ToString() + ".mol";

            lock (noneparallel)
                conv.WriteFile(mol, Path.Combine("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\examples", filename));
            string minimizeOutput;

            using (Process proc = new Process())
            {
                // C:\Program Files (x86)\OpenBabel-3.1.1

                proc.StartInfo.FileName = "C:\\Program Files (x86)\\OpenBabel-3.1.1\\obminimize.exe";
                proc.StartInfo.Arguments = "-ff UFF " + filename;

                //proc.StartInfo.Arguments = "-n200 minimize.mol"; //can add arguments here like number of iterations,
                // or '-c' convergence criteria
                proc.StartInfo.ErrorDialog = false;
                proc.StartInfo.WorkingDirectory = "C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\examples";
                //proc.StartInfo.RedirectStandardError = true;
                //proc.StartInfo.UseShellExecute = false;
                proc.StartInfo.RedirectStandardOutput = true;
                //proc.StartInfo.RedirectStandardInput = false;

                Console.Write("starting OBMinimize...");
                proc.Start();

                minimizeOutput = proc.StandardOutput.ReadToEnd();
                proc.WaitForExit();

            }
            lock (noneparallel)
                conv.ReadString(mol, minimizeOutput);

            //File.Delete("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\examples\\"+ filename);




            return mol;

        }






    }

}






