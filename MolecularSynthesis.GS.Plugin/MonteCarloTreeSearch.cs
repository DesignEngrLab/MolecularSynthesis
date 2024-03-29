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

namespace MolecularSynthesis.GS.Plugin
{
    public class MCTS : SearchProcess
    {
        // give desiredMoment
        // [] desiredMoment = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        static double[] desiredLenghtAndRadius = new double[] { 300, 50 };
        static Random rnd = new Random(0);

        public MCTS(GlobalSettings settings) : base(settings)
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
            get { return "MCTS"; }
        }
        protected override void Run()
        {
            //var candidates = new SimplePriorityQueue<candidate, double>();

            // generate a random number 0 or 1 to decide the next rule is from RS0 or RS1
            //Random rnd = new Random();
            //rnd.Next(0, 2); // generate 0 or 1

            // use 10000 is that DS use 3000-70000 iteration for 9*9 go play , so guess 10000 is enough
            int iteration = 10000;
            //TreeCandidate node1 = new TreeCandidate() { S = 0, n=0, UCB=0 };

            // 1. check if this is the leaf node, if not, go to step 2 until it is a leaf node,if yes go to step 3
            // 2. find the children who has the best UCB value
            // 3. do random simulation
            // 4. update S,n,UCB value for the whole tree

<<<<<<< HEAD
            //when the search is within range of <> , BFS is better
            //when the search is within range of <> , MCTS is better

            // rollout process should match linker length , add a collector to collect solutions from trees with different depth
            // change for loop to while loop, set stop criteria, like n=50
            // careful for result from evaluation, should include posive value and negative value
            
=======
>>>>>>> parent of 726848f (Still testing MCTS)
            TreeCandidate StartState = new TreeCandidate(seedCandidate);

            StartState.S = 0;
            StartState.n = 0;
            StartState.UCB = double.MaxValue;
            StartState.Children = new List<TreeCandidate>();
            StartState.Parent = null;

            int IterationTimes = 0;


            // while (current.n<50)
            for (int i = 0; i < iteration; i++)
            {
                // if abs(current.S - target value)  < stop criteria 
                //  record this recipe

                // need to save S value and n value, delete the added graph, back to StartState                                                  
                TreeCandidate current = StartState;
                while (current.Children.Count > 0)
                    current = SelectPromisingNode(current);// until at leaf node               

                if (current.n == 0)
                    current.S = Rollout(current);
                else
                {
                    // add all possible actions under one parent node
<<<<<<< HEAD
                    int RS0 = 0;
                    foreach (var option in current.recipe)
                    {
                        // need to recognize how many Rules from Ruleset0 exist
                        if (option.ruleSetIndex == 0)
                            RS0++;
                    }
                    // go RS0 RS2 RS1 
                    if (RS0 < 5)
=======
                    if (current.recipe.Count < 3)
>>>>>>> parent of 726848f (Still testing MCTS)
                    {
                        AddNewNode(current);
                        SearchIO.output("Children number = " + current.Children.Count + "**********");
                        current = SelectPromisingNode(current);
                        current.S = Rollout(current);
                    }
                    else
                        current.S = Rollout(current);
                }
                BackPropogation(FindAllParents(current), current);
                IterationTimes = DisplayData(IterationTimes, MCTSProcess, current);
            }
            ReportFinalData(StartState, MCTSProcess);

<<<<<<< HEAD
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
=======
                IterationTimes = IterationTimes + 1;

                SearchIO.output("Iteration times: " + IterationTimes);
                SearchIO.output("S = " + current.S);
                SearchIO.output("n = " + current.n);
                SearchIO.output("Children number = " + current.Children.Count + "*************");
                SearchIO.output("number of recipe: " + current.recipe.Count);
                SearchIO.output("current node recipe:");

                if (current.recipe.Count == 0)
                {
                    SearchIO.output("no recipe" + "------------");
                }
                else
                {
                    foreach (var option in current.recipe)
                    {
                        //SearchIO.output(current.recipe[j].ruleSetIndex + " " + current.recipe[j].optionNumber);
                        SearchIO.output(option.ruleSetIndex + " " + option.ruleNumber + "------------");

                    }
                }

                SearchIO.output("-------------------------------------------------------------------------");
            }

            //TreeCandidate seed = new TreeCandidate(seedCandidate);
            var FinalResult = FinalRecipe(StartState);
            SearchIO.output("Solution is down below: ");
            foreach (var option in FinalResult.recipe)
            {
                SearchIO.output(option.ruleSetIndex + " " + option.ruleNumber + "------------");
            }


>>>>>>> parent of 726848f (Still testing MCTS)

            string seperateline = "-------------------------------------------------------------------------";
            MCTSProcess.Add(seperateline);
            return IterationTimes;
        }

        public double CalculateUcb(TreeCandidate child)
        {
            if (child.n == 0)
                return double.MaxValue;
            else
                return child.S / child.n + 200 * Math.Sqrt(Math.Log(child.Parent.n) / child.n);
        }

            //create the bestchild as an intermidiate variable
        public TreeCandidate SelectPromisingNode(TreeCandidate current)
        {
            TreeCandidate bestChild = null;
<<<<<<< HEAD
=======


>>>>>>> parent of 726848f (Still testing MCTS)
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
                    
                    if (child.S > bestS)
                    {
                        bestS = child.S;
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

        public double Rollout(TreeCandidate candidate)
        {
            double score;
            int RS0 = candidate.ruleSetIndicesInRecipe.Count(rsIndex => rsIndex == 0);
            int RS1 = 0;
            TreeCandidate child = (TreeCandidate)candidate.copy();                     

<<<<<<< HEAD
            // (RS0 R1) 
            while (RS0 < 5)
=======
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

                //var options = rulesets[0].recognize(current.graph)
            }

            while (RS0 < 3)
>>>>>>> parent of 726848f (Still testing MCTS)
            {
                var option0 = rulesets[0].recognize(child.graph);
                int WhichRuleset = rnd.Next(0, 2);
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
                }
            }
            var option2 = rulesets[2].recognize(child.graph);
<<<<<<< HEAD
=======
            
>>>>>>> parent of 726848f (Still testing MCTS)
            if (option2.Count != 1)
                Console.WriteLine("how?!?"+"|||  Option2="+ option2.Count);
            option2[0].apply(child.graph, null);
            child.addToRecipe(option2[0]);

<<<<<<< HEAD
            return (double)Evaluation.TotalAtomMass(child);
=======
            //if (!candidate.recipe.Contains(option2[0]))
            //{
            //    option2[0].apply(candidate.graph, null);
            //}

            // use openbabel for evaluation
            OBMol resultMol = OBFunctions.designgraphtomol(child.graph);
            resultMol = OBFunctions.InterStepMinimize(resultMol);
            OBFunctions.updatepositions(child.graph, resultMol);

            score = -Evaluation.distance(child, desiredLenghtAndRadius);
            return score;

            // just for testing , no need for openbabel
            //double TotalMass = Evaluation.TotalAtomMass(child);
            //return TotalMass;



>>>>>>> parent of 726848f (Still testing MCTS)
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






