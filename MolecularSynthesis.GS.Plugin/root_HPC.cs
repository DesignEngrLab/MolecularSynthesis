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
using System.Threading.Tasks;

namespace MolecularSynthesis.GS.Plugin
{
    public class MCTS_Parallel_root_HPC : SearchProcess
    {
        // give desiredMoment
        // [] desiredMoment = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        // RS0 R4 R5 R6 R7; RS1 R8 R9 target
        // RS0 R4 R2 R4 R6 R6  solution 897.486 283.802
        // RS0 R4 R5 R6 R2 R3 894.260 288.637
        static double[] desiredLenghtAndRadius = new double[] { 891.222, 288.221 };
        static Random rnd = new Random(0);
        // the linker without any rules will have 41.675 in length and 90.584 in radius 
        // 891.222-41.675=849.5
        // 288.221- 90.584 = 197.635
        // 872.18


        public MCTS_Parallel_root_HPC(GlobalSettings settings) : base(settings)
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
            get { return "MCTS_root_HPC"; }
        }


        protected override void Run()
        {

            // add a timer to compare running time for different paralellization method
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();


            //var candidates = new SimplePriorityQueue<candidate, double>();

            // generate a random number 0 or 1 to decide the next rule is from RS0 or RS1
            //Random rnd = new Random();
            //rnd.Next(0, 2); // generate 0 or 1

            // use 10000 is that DS use 3000-70000 iteration for 9*9 go play , so guess 10000 is enough
            int iteration = 10;


            //TreeCandidate node1 = new TreeCandidate() { S = 0, n=0, UCB=0 };

            // 1. check if this is the leaf node, if no go to step 2 until it is a leaf node,if yes go to step 3
            // 2. find the children who has the best UCB value
            // 3. do random simulation
            // 4. update S,n,UCB value for the whole tree

            //when the search is within range of <> , BFS is better
            //when the search is within range of <> , MCTS is better

            var stopwatch = new Stopwatch();

            TreeCandidate StartState = new TreeCandidate(seedCandidate);

            StartState.S = 0;
            StartState.n = 0;
            StartState.UCB = double.MaxValue;
            StartState.Children = new List<TreeCandidate>();
            StartState.Parent = null;

            int IterationTimes = 0;
            List<string> MCTSProcess = new List<string>();

            var current = StartState;

            // 1. leaf parallel 2. node parallel 3 tree parallel 

            //Parallel.ForEach(list_lines, line =>
            //{
            //    //Your stuff
            //});


            //while (current.n < 10)           
            for (int i = 0; i <= iteration; i++)
            {
                // need to save S value and n value, delete the added graph, back to StartState                                                  

                current = StartState;

                while (current.Children.Count > 0)
                {
                    current = SelectPromisingNode(current);// until at leaf node               

                }

                if (current.n == 0)
                {
                    current.S = Rollout(current);
                }


                else
                {
                    // add all possible actions under one parent node
                    int RS0 = 0;
                    foreach (var option in current.recipe)
                    {
                        // need to recognize how many Rules from Ruleset0 exist
                        if (option.ruleSetIndex == 0)
                        {
                            RS0 = RS0 + 1;
                        }
                    }
                    // go RS0 RS2 RS1 
                    if (RS0 < 5)
                    {
                        AddNewNode(current);
                        string ChildrenInformation = "Children number = " + current.Children.Count.ToString() + "**********";
                        MCTSProcess.Add(ChildrenInformation);
                        current = SelectPromisingNode(current);
                        current.S = Rollout(current);
                    }
                    else
                    {
                        current.S = Rollout(current);
                    }

                }



                BackPropogation(FindAllParents(current), current);

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
            };

            //TreeCandidate seed = new TreeCandidate(seedCandidate);
            var FinalResult = FinalRecipe(StartState);
            string SolutionIsDownBelow = "Solution is down below: ";
            MCTSProcess.Add(SolutionIsDownBelow);

            foreach (var option in FinalResult.recipe)
            {
                string SolutionInformation = option.ruleSetIndex + " " + option.ruleNumber + "------------";
                MCTSProcess.Add(SolutionInformation);
            }

            stopWatch.Stop();
            // Get the elapsed time as a TimeSpan value.
            TimeSpan ts = stopWatch.Elapsed;
            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
            ts.Hours, ts.Minutes, ts.Seconds,
            ts.Milliseconds / 10);

            MCTSProcess.Add(elapsedTime);

            System.IO.File.WriteAllLines(@"C:\Users\zhang\source\repos\MolecularSynthesis\output\MCTSProcessRecord.txt", MCTSProcess);

            
        }

        public double CalculateUcb(TreeCandidate child)
        {
            if (child.n == 0)
                return double.MaxValue;
            else
                return child.S / child.n + 500 * Math.Sqrt(Math.Log(child.Parent.n) / child.n);
        }

        public TreeCandidate SelectPromisingNode(TreeCandidate current)
        {
            //create the bestchild as an intermidiate variable

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

            //return SelectPromisingNode(bestChild);
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
            current.n = current.n + 1;

            foreach (TreeCandidate treeCandidate in parentpath)
            {

                treeCandidate.n++;
                treeCandidate.S += current.S;

            }
        }

        public double Rollout(TreeCandidate candidate)
        {
            double score;
            int RS0 = 0;
            int RS1 = 0;
            TreeCandidate child = (TreeCandidate)candidate.copy();

            //option
            //public int ruleSetIndex { get; set; }
            //public int optionNumber { get; set; }
            foreach (var option in child.recipe)
            {
                // need to recognize how many Rules from Ruleset0 exist
                if (option.ruleSetIndex == 0)
                {
                    RS0 = RS0 + 1;
                }
            }

            while (RS0 < 5)
            {
                //rnd.Next(0, 2); // generate 0 or 1

                var option0 = rulesets[0].recognize(child.graph);
                int WhichRuleset = rnd.Next(0, 2);
                //int WhichRuleset = 0;

                if (WhichRuleset == 0)
                {
                    RS0 = RS0 + 1;
                    if (option0.Count > 0)
                    {
                        option0 = rulesets[0].recognize(child.graph);
                        var Randomoption0 = rnd.Next(0, option0.Count);
                        option0[Randomoption0].apply(child.graph, null);
                        // dont need to add options into recipe in rollout process
                        child.addToRecipe(option0[Randomoption0]);
                    }
                }

                else
                {
                    var option1 = rulesets[1].recognize(child.graph);
                    if (option1.Count > 0)
                    {
                        RS1 = RS1 + 1;
                        var Randomoption1 = rnd.Next(0, option1.Count);
                        option1[Randomoption1].apply(child.graph, null);
                        child.addToRecipe(option1[Randomoption1]);
                    }

                    //candidate=RecognizeChooseApply.GenerateAllNeighbors(current, rulesets, false, false, true)

                }
            }

            var option2 = rulesets[2].recognize(child.graph);

            if (option2.Count != 1)
                Console.WriteLine("how?!?" + "|||  Option2=" + option2.Count);
            option2[0].apply(child.graph, null);
            child.addToRecipe(option2[0]);

            //if (!candidate.recipe.Contains(option2[0]))
            //{
            //    option2[0].apply(candidate.graph, null);
            //}

            // use openbabel for evaluation
            //OBMol resultMol = OBFunctions.designgraphtomol(child.graph);
            //resultMol = OBFunctions.InterStepMinimize(resultMol);
            //OBFunctions.updatepositions(child.graph, resultMol);

            //score = -Evaluation.distance(child, desiredLenghtAndRadius);
            //return score;

            // just for testing , no need for openbabel
            double TotalMass = -Evaluation.TotalAtomMass(child);
            return TotalMass;



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
    }

}






