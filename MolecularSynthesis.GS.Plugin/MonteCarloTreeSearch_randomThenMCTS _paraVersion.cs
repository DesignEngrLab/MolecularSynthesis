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
using System.Diagnostics;
using System.Threading;
using System.Threading.Tasks;

namespace MolecularSynthesis.GS.Plugin
{
    public class MonteCarloTreeSearch_randomThenMCTS_paraVersion : SearchProcess
    {
        // give desiredMoment
        // [] desiredMoment = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        // RS0 R7 R3; RS1 R1 R2
        //RS0 R3 L=134;
        //static double[] desiredLenghtAndRadius = new double[] { 245.277, 89.53 };
        static Random rnd = new Random();

        static double[] desiredLenghtAndRadius = new double[] { 565, 140 };
        // RS0 R3 R4 R5 R6
        static TreeCandidate noneparallel = new TreeCandidate(new candidate());

        public MonteCarloTreeSearch_randomThenMCTS_paraVersion(GlobalSettings settings) : base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 2;
            AutoPlay = true;
        }

        /// <summary>
        ///   Gets the text describing that is displayed in the menu. It must be overridden in the methods.
        /// </summary>
        /// <value>The text.</value>
        public override string text
        {
            get { return "MonteCarloTreeSearch_randomThenMCTS_paraVersion"; }
        }
        protected override void Run()
        {
            //var candidates = new SimplePriorityQueue<candidate, double>();

            // generate a random number 0 or 1 to decide the next rule is from RS0 or RS1
            //Random rnd = new Random();
            //rnd.Next(0, 2); // generate 0 or 1

            // use 10000 is that DS use 3000-70000 iteration for 9*9 go play , so guess 10000 is enough
            int iteration = 10;

            //TreeCandidate node1 = new TreeCandidate() { S = 0, n=0, UCB=0 };

            // 1. check if this is the leaf node, if not, go to step 2 until it is a leaf node,if yes go to step 3
            // 2. find the children who has the best UCB value
            // 3. do random simulation
            // 4. update S,n,UCB value for the whole tree


            // rollout process should match linker length , add a collector to collect solutions from trees with different depth
            // change for loop to while loop, set stop criteria, like n=50
            // careful for result from evaluation, should include posive value and negative value

            var timer = new Stopwatch();
            timer.Start();

            Parallel.For(0, 10, count =>
            {


                TreeCandidate StartState = new TreeCandidate(seedCandidate);

                StartState.S = 0;
                StartState.n = 0;
                StartState.UCB = double.MaxValue;
                StartState.Children = new List<TreeCandidate>();
                StartState.Parent = null;

                int IterationTimes = 0;
                List<string> MCTSProcess = new List<string>();
                List<string> resultCollector = new List<string>();

                double[] Everystep = new double[2];
                double score = 0;



                var rand = new Random();

                //TreeCandidate current = StartState;

                //while (current.n<50)
                for (int i = 0; i < iteration; i++)
                {
                    Console.WriteLine("-----------------------------iterationtime=" + i);
                    // if abs(current.S - target value)  < stop criteria 
                    //  record this recipe

                    // need to save S value and n value, delete the added graph, back to StartState                                                  
                    //TreeCandidate current = (TreeCandidate)StartState.copy();
                    TreeCandidate current = StartState;

                    // Random search


                    //------------------------------------------------------------------
                    // for random search here, first we need to determine the level of candidate in the tree, then pick a random leaf at each level  

                    if (i < 1000)
                    {
                        var Levelnumber = rand.Next(1, 15);

                        var Leafnumber = 0;
                        var n = 0;
                        while (n < Levelnumber)
                        {
                            if (current.Children.Count == 0)
                            {
                                AddNewNode(current);

                            }
                            else
                            {
                                Leafnumber = rand.Next(0, current.Children.Count);
                                current = (TreeCandidate)current.Children[Leafnumber].copy();
                                n = n + 1;
                            }
                            //var resultMoltest = OBFunctions.designgraphtomol(current.graph);
                            //resultMoltest = justMinimize(resultMoltest);
                        }


                        var option2 = rulesets[2].recognize(current.graph);
                        option2[0].apply(current.graph, null);
                        current.addToRecipe(option2[0]);

                        var resultMol = OBFunctions.designgraphtomol(current.graph);
                        resultMol = justMinimize(resultMol);
                        OBFunctions.updatepositions(current.graph, resultMol);

                        score = Evaluation.distance(current, desiredLenghtAndRadius);
                        BackPropogation(FindAllParents(current), current);
                    }

                    //--------------------------------------------------------------------

                    else
                    {
                        // MCTS
                        while (current.Children.Count > 0)
                            current = SelectPromisingNode(current);// until at leaf node               

                        if (current.n == 0)
                        {
                            Everystep = Rollout(current);
                            current.S = Everystep[0];
                            score = Everystep[1];
                        }
                        else
                        {
                            // add all possible actions under one parent node
                            int RS0 = 0;
                            foreach (var option in current.recipe)
                            {
                                // need to recognize how many Rules from Ruleset0 exist
                                if (option.ruleSetIndex == 0)
                                    RS0++;
                            }
                            // go RS0 RS2 RS1 
                            if (RS0 < 5)
                            {
                                AddNewNode(current);
                                string ChildrenInformation = "Children number = " + current.Children.Count.ToString() + "**********";
                                MCTSProcess.Add(ChildrenInformation);
                                current = SelectPromisingNode(current);


                                Everystep = Rollout(current);
                                current.S = Everystep[0];
                                score = Everystep[1];

                            }
                            else
                            {
                                Everystep = Rollout(current);
                                current.S = Everystep[0];
                                score = Everystep[1];
                            }
                        }
                    }
                    // --------------------collect current evaluation value at each iteration-----------------
                    //var resultMol = OBFunctions.designgraphtomol(current.graph);
                    //resultMol = justMinimize(resultMol);
                    //OBFunctions.updatepositions(current.graph, resultMol);

                    //var score = Evaluation.distance(current, desiredLenghtAndRadius);
                    resultCollector.Add(score.ToString());

                    //--------------------------------------------------------------------------------------

                    BackPropogation(FindAllParents(current), current);
                    //IterationTimes = DisplayData(IterationTimes, MCTSProcess, current);
                }

                //ReportFinalData(StartState, MCTSProcess);

                var filename = "MCTSthennRandom_data_foruse" + Thread.CurrentThread.ManagedThreadId.ToString();
                filename = filename + ".txt";
                System.IO.File.WriteAllLines(@"C:\Users\zhang\source\repos\MolecularSynthesis\examples\" + filename, resultCollector);
            });

            timer.Stop();
            TimeSpan ts = timer.Elapsed;

            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
            ts.Hours, ts.Minutes, ts.Seconds,
            ts.Milliseconds / 10);
            Console.WriteLine("RunTime " + elapsedTime);



        }

        private void ReportFinalData(TreeCandidate StartState, List<string> MCTSProcess)
        {


            //TreeCandidate seed = new TreeCandidate(seedCandidate);
            var FinalResult = FinalRecipe(StartState);
            string SolutionIsDownBelow = "Solution is down below: ";
            MCTSProcess.Add(SolutionIsDownBelow);

            foreach (var option in FinalResult.recipe)
            {
                string SolutionInformation = option.ruleSetIndex + " " + option.ruleNumber + "------------";
                MCTSProcess.Add(SolutionInformation);
            }

            System.IO.File.WriteAllLines(@"C:\Users\zhang\source\repos\MolecularSynthesis\output\MCTSProcessRecord.txt", MCTSProcess);
        }

        private int DisplayData(int IterationTimes, List<string> MCTSProcess, TreeCandidate current)
        {
            IterationTimes = IterationTimes + 1;
            SearchIO.output("IterationTimes = ", IterationTimes);

            string times = "Iteration times: " + IterationTimes.ToString();
            MCTSProcess.Add(times);
            string CurrentSValue = "S = " + current.S.ToString();
            MCTSProcess.Add(CurrentSValue);
            string CurrentnValue = "n = " + current.n.ToString();
            MCTSProcess.Add(CurrentnValue);

            if (IterationTimes > 1)
            {

                string CurrentUCBValue = "UCB = " + CalculateUcb(current).ToString();
                MCTSProcess.Add(CurrentUCBValue);
            }
            string ChildrenNumber = "Children number = " + current.Children.Count.ToString() + "*************";
            MCTSProcess.Add(ChildrenNumber);
            string NumberOfRecipe = "number of recipe: " + current.recipe.Count.ToString();
            MCTSProcess.Add(NumberOfRecipe);
            string CurrentNodeRecipe = "current node recipe:";
            MCTSProcess.Add(CurrentNodeRecipe);

            if (current.recipe.Count == 0)
            {
                string NoRecipe = "no recipe" + "------------";
                MCTSProcess.Add(NoRecipe);
            }
            else
            {
                foreach (var option in current.recipe)
                {
                    //SearchIO.output(current.recipe[j].ruleSetIndex + " " + current.recipe[j].optionNumber);
                    string OptionInformation = option.ruleSetIndex.ToString() + " " + option.ruleNumber.ToString() + "------------";
                    MCTSProcess.Add(OptionInformation);
                }
            }

            string seperateline = "-------------------------------------------------------------------------";
            MCTSProcess.Add(seperateline);
            return IterationTimes;
        }

        public double CalculateUcb(TreeCandidate child)
        {
            if (child.n == 0)
                return double.MaxValue;
            else
                return child.S / child.n + 20 * Math.Sqrt(Math.Log(child.Parent.n) / child.n);
        }

        //create the bestchild as an intermidiate variable
        public TreeCandidate SelectPromisingNode(TreeCandidate current)
        {
            TreeCandidate bestChild = null;
            while (current.Children.Count != 0)
            {
                double bestUcb = double.MinValue;
                foreach (TreeCandidate child in current.Children)
                {
                    double Ucb = CalculateUcb(child);
                    child.UCB = Ucb;
                    if (Ucb > bestUcb)
                    {
                        bestUcb = Ucb;
                        bestChild = child;
                    }
                }
                current = bestChild;
            }
            return bestChild;
        }

        public TreeCandidate FinalRecipe(TreeCandidate StartState)
        {
            TreeCandidate bestChild = null;
            while (StartState.Children.Count != 0)
            {
                double bestS = double.MinValue;
                foreach (TreeCandidate child in StartState.Children)
                {

                    if (child.S / child.n > bestS)
                    {
                        bestS = child.S / child.n;
                        bestChild = child;
                    }
                }
                StartState = bestChild;

            }

            return bestChild;
        }


        public void AddNewNode(TreeCandidate current)
        {
            // need to add one avaiable option from current ,add options into recipe

            var option0 = rulesets[0].recognize(current.graph);
            //option0.AddRange(rulesets[1].recognize(current.graph));
            var option1 = rulesets[1].recognize(current.graph);
            int PotenialOptionNumber = option1.Count + option0.Count;
            //int PotenialOptionNumber = option0.Count;


            //var option0 = rulesets[0].recognize(candidate.graph);
            for (int i = 0; i < PotenialOptionNumber; i++)
            {
                var child = (TreeCandidate)current.copy();
                child.Parent = current;
                child.Children = new List<TreeCandidate>();
                child.n = 0;
                child.S = 0;
                child.UCB = double.MaxValue;

                if (i < option0.Count)
                {
                    option0 = rulesets[0].recognize(child.graph);
                    option0[i].apply(child.graph, null);
                    child.addToRecipe(option0[i]);
                }
                else
                {
                    option1 = rulesets[1].recognize(child.graph);
                    option1[i - option0.Count].apply(child.graph, null);
                    child.addToRecipe(option1[i - option0.Count]);
                }

                current.Children.Add(child);


            }




        }

        public void BackPropogation(List<TreeCandidate> parentpath, TreeCandidate current)
        {
            current.n++;
            foreach (TreeCandidate treeCandidate in parentpath)
            {
                treeCandidate.n++;
                treeCandidate.S += current.S;
            }
        }

        public double[] Rollout(TreeCandidate candidate)
        {
            double score;
            int RS0 = candidate.ruleSetIndicesInRecipe.Count(rsIndex => rsIndex == 0);
            int RS1 = 0;
            var rand = new Random();
            TreeCandidate child = (TreeCandidate)candidate.copy();



            //-----------------------new rollout method-------------------------------------------
            // as long as RS0<=5 , RS0 can be 1 2 3 4 5, number of RS1 is random

            var numberOFBB = rnd.Next(RS0 + 1, 6);

            // generate backbone group rules
            while (RS0 < numberOFBB)

            {
                RS0 = RS0 + 1;
                var option0 = rulesets[0].recognize(child.graph);
                var totalOptionNumber = option0.Count();
                var OptionNumber = rnd.Next(0, totalOptionNumber);
                option0[OptionNumber].apply(child.graph, null);

            }

            var option1 = rulesets[1].recognize(child.graph);
            var option2 = rulesets[2].recognize(child.graph);

            // generate functional group rules
            if (option1.Count == 0)
            {
                option2 = rulesets[2].recognize(child.graph);
                option2[0].apply(child.graph, null);
            }

            else
            {
                var PossibleFGnumber = option1.Count;

                // check how many position are available
                var totalFGnumber = rand.Next(0, PossibleFGnumber / 9);


                for (int j = 0; j < totalFGnumber; j++)
                {
                    var OptionNumber2 = rand.Next(0, PossibleFGnumber);


                    option1[OptionNumber2].apply(child.graph, null);
                    child.addToRecipe(option1[OptionNumber2]);

                    option1 = rulesets[1].recognize(child.graph);
                    PossibleFGnumber = option1.Count;
                }

                option2 = rulesets[2].recognize(child.graph);
                option2[0].apply(child.graph, null);
                //candidate.addToRecipe(option2[0]);

            }


            //-----------------------------old rollout method----------------------------------------
            //var lengthofbb = rnd.next(rs0,6);
            //// (rs0 r1) 
            //while (rs0 < lengthofbb)
            //{
            //    var option0 = rulesets[0].recognize(child.graph);
            //    int whichruleset = rnd.next(0, 2);
            //    if (whichruleset == 0)
            //    {
            //        rs0 = rs0 + 1;
            //        if (option0.count > 0)
            //        {
            //            option0 = rulesets[0].recognize(child.graph);
            //            var randomoption0 = rnd.next(0, option0.count);
            //            option0[randomoption0].apply(child.graph, null);
            //            // dont need to add options into recipe in rollout process
            //            child.addtorecipe(option0[randomoption0]);
            //        }
            //    }
            //    else
            //    {
            //        var option1 = rulesets[1].recognize(child.graph);
            //        if (option1.count > 0)
            //        {
            //            rs1 = rs1 + 1;
            //            var randomoption1 = rnd.next(0, option1.count);
            //            option1[randomoption1].apply(child.graph, null);
            //            child.addtorecipe(option1[randomoption1]);
            //        }
            //    }
            //}


            //-----------------stop Growing--------------------------------------------------
            //var option2 = rulesets[2].recognize(child.graph);
            //if (option2.Count != 1)
            //    Console.WriteLine("how?!?" + "|||  Option2=" + option2.Count);
            //option2[0].apply(child.graph, null);
            //child.addToRecipe(option2[0]);

            OBMol resultMol = OBFunctions.designgraphtomol(child.graph);
            resultMol = justMinimize(resultMol);
            OBFunctions.updatepositions(child.graph, resultMol);

            score = Evaluation.distance(child, desiredLenghtAndRadius);

            double[] ResultInformation = new double[2];
            ResultInformation[0] = 100000 - score;
            ResultInformation[1] = score;

            return ResultInformation;




            //return (double)Evaluation.TotalAtomMass(child);
        }
        public List<TreeCandidate> FindAllParents(TreeCandidate current)
        {
            var parents = new List<TreeCandidate>();
            TreeCandidate child = (TreeCandidate)current.copy();
            //child.Parent = current;
            //.Parent = current.Parent;
            //TreeCandidate child = current;

            //int height = GetHeight(child);

            while (child.Parent != null)
            //for (int i = 0; i < height; i++)
            {
                parents.Add(child.Parent);
                child = child.Parent;
            }

            return parents;
        }

        public int GetHeight(TreeCandidate current)
        {
            int height = 0;


            TreeCandidate child = (TreeCandidate)current.copy();
            //child.Parent = current;
            //child.Parent = current.Parent;

            while (child.Parent != null)
            {
                height = height + 1;
                child = child.Parent;
            }

            return height;
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






