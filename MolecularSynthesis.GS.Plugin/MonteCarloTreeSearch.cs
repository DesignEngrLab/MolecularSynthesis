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

        //public class TreeNode<T>
        //{
        //    public T Data { get; set; }
        //    public TreeNode<T> Parent { get; set; }

        //    public List<TreeNode<T>> Children { get; set; }
        //    public int GetHeight()
        //    {
        //        int height = 1;
        //        TreeNode<T> current = this;
        //        while (current.Parent != null)
        //        {
        //            height++;
        //            current = current.Parent;

        //        }
        //        return height;
        //    }
        //}

        //public class Tree<T>
        //{            
        //    public TreeNode<T> Root { get; set; }
        //}

        //class program
        //{
        //    static void Main(string[] args)
        //    {
        //        Tree<int> tree = new Tree<int>();
        //        tree.Root = new TreeNode<int>() { Data = 100 };
        //        tree.Root.Children = new List<TreeNode<int>>
        //        {
        //            new TreeNode<int>(){ Data=50, Parent=tree.Root },
        //            new TreeNode<int>(){ Data=1, Parent=tree.Root },
        //            new TreeNode<int>(){ Data=150, Parent=tree.Root}
        //        };
        //        tree.Root.Children[2].Children = new List<TreeNode<int>>()
        //        {
        //            new TreeNode<int>()
        //            { Data=30, Parent= tree.Root.Children[2]}
        //        };
        //    }

        //}

        protected override void Run()
        {
            //var candidates = new SimplePriorityQueue<candidate, double>();

            // generate a random number 0 or 1 to decide the next rule is from RS0 or RS1
            //Random rnd = new Random();
            //rnd.Next(0, 2); // generate 0 or 1

            // use 10000 is that DS use 3000-70000 iteration for 9*9 go play , so guess 10000 is enough
            int iteration = 2000;
            //TreeCandidate node1 = new TreeCandidate() { S = 0, n=0, UCB=0 };

            // 1. check if this is the leaf node, if no go to step 2 until it is a leaf node,if yes go to step 3
            // 2. find the children who has the best UCB value
            // 3. do random simulation
            // 4. update S,n,UCB value for the whole tree

            TreeCandidate StartState = new TreeCandidate(seedCandidate);

            StartState.S = 0;
            StartState.n = 0;
            StartState.UCB = double.MaxValue;
            StartState.Children = new List<TreeCandidate>();
            StartState.Parent = null;

            int IterationTimes = 0;


            var current = StartState;

            for (int i = 0; i < iteration; i++)
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
                    if (current.recipe.Count < 3)
                    {
                        AddNewNode(current);
                        SearchIO.output("Children number = " + current.Children.Count + "**********");
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



        }
        public double CalculateUcb(TreeCandidate child)
        {
            if (child.n == 0)
                return double.MaxValue;
            else
                return child.S / child.n + 200 * Math.Sqrt(Math.Log(child.Parent.n) / child.n);
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
            // option1 = rulesets[1].recognize(current.graph);
            //int PotenialOptionNumber = option1.Count + option0.Count;
            int PotenialOptionNumber = option0.Count;


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
                //else
                //{
                    //option1 = rulesets[1].recognize(child.graph);
                //    option1[i - option0.Count].apply(child.graph, null);
                //    child.addToRecipe(option1[i - option0.Count]);
                //}

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

                //var options = rulesets[0].recognize(current.graph)
            }

            while (RS0 < 3)
            {
                //rnd.Next(0, 2); // generate 0 or 1

                var option0 = rulesets[0].recognize(child.graph);
                //int WhichRuleset = rnd.Next(0, 2);
                int WhichRuleset = 0;

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
                Console.WriteLine("how?!?"+"|||  Option2="+ option2.Count);
            option2[0].apply(child.graph, null);
            child.addToRecipe(option2[0]);

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






