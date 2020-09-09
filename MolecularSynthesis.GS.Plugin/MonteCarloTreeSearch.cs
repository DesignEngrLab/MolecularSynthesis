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
            int iteration = 10000;
            //TreeCandidate node1 = new TreeCandidate() { S = 0, n=0, UCB=0 };


            // 1. check if this is the leaf node, if no go to step 2 until it is a leaf node,if yes go to step 3
            // 2. find the children who has the best UCB value
            // 3. do random simulation
            // 4. update S,n,UCB value for the whole tree
            TreeCandidate current = (TreeCandidate)seedCandidate;
            current.S = 0;
            current.n = 0;
            current.UCB = double.MaxValue;
            current.Children = new List<TreeCandidate>();
            //current.option = 0;

            for (int i = 0; i < iteration; i++)
            {
                if (current.Children == null)
                {
                    if (current.n == 0)
                    {
                        Rollout(current);
                    }
                    else
                    {
                        current = AddNewNode(current);
                        Rollout(current);

                    }
                }
                else
                {
                    current = SelectPromisingNode(current);
                }

                //current = SelectPromisingNode(current);

                //Rollout(current);

                BackPropogation(FindAllParents(current));

            }

            TreeCandidate seed = (TreeCandidate)seedCandidate;
            var FinalResult=SelectPromisingNode(seed);
            SearchIO.output(FinalResult.recipe);
        }
        public double CalculateUcb(TreeCandidate child)
        {
            return child.S / (double)child.n + 2 * Math.Sqrt(Math.Log(child.Parent.n) / child.n);
        }
        public TreeCandidate SelectPromisingNode(TreeCandidate current)
        {
            TreeCandidate bestChild = null;
            double bestUcb = double.MinValue;

            while (current.Children != null)
            {
                foreach (TreeCandidate child in current.Children)
                {
                    double Ucb = CalculateUcb(child);
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

        public TreeCandidate AddNewNode(TreeCandidate current)
        {
            //is this right?
            var child = new TreeCandidate(current);
            child.Parent = current;
            child.Children = null;
            child.n = 0;
            child.S = 0;
            child.UCB = double.MinValue;

            current.Children.Add(child);
            return child;
        }

        public void BackPropogation(List<TreeCandidate> parentpath)
        {
            foreach (var treeCandidate in parentpath)
            {
                treeCandidate.n++;
                treeCandidate.S++;

            }
        }

        public double Rollout(TreeCandidate candidate)
        {
            double score;
            int RS0 = 0;
            var options = rulesets;
            //candidate.ruleSetIndicesInRecipe();
            //var childrenCandidate = RecognizeChooseApply.GenerateAllNeighbors(current, rulesets, false, false, true);
            //var a = candidate.recipe[1];

            //option
            //public int ruleSetIndex { get; set; }
            //public int optionNumber { get; set; }
            foreach (var option in candidate.recipe)
            {
                // need to recognize how many Rules from Ruleset0 exist
                if (option.ruleSetIndex == 0)
                {
                    RS0 = RS0 + 1;
                }

                //var options = rulesets[0].recognize(current.graph)
            }

            var option0 = rulesets[0].recognize(candidate.graph);
            var option1 = rulesets[1].recognize(candidate.graph);
            var option2 = rulesets[2].recognize(candidate.graph);

            while (RS0 != 5)
            {
                Random rnd = new Random();
                //rnd.Next(0, 2); // generate 0 or 1
                if (rnd.Next(0, 2) == 0)
                {
                    RS0 = RS0 + 1;
                    option0[rnd.Next(option0.Count)].apply(candidate.graph, null);
                    candidate.addToRecipe(option0[option0.Count]);
                }
                else
                {
                    option1[rnd.Next(option1.Count)].apply(candidate.graph, null);
                    candidate.addToRecipe(option1[option1.Count]);
                }

                //candidate=RecognizeChooseApply.GenerateAllNeighbors(current, rulesets, false, false, true)
            }
            option2[0].apply(candidate.graph, null);
            var resultMol = OBFunctions.designgraphtomol(candidate.graph);
            resultMol = OBFunctions.InterStepMinimize(resultMol);
            OBFunctions.updatepositions(seedGraph, resultMol);

            // distance is less the better, do we need to add a negative before score 
            score = -Evaluation.distance(candidate, desiredLenghtAndRadius);
            return score;
        }

        public List<TreeCandidate> FindAllParents(TreeCandidate current)
        {
            var parents = new List<TreeCandidate>();
            int height = GetHeight(current);

            for (int i = 0; i < height; i++)
            {
                parents.Add(current.Parent);
            }

            return parents;
        }

        public int GetHeight(TreeCandidate current)
        {
            int height = 0;
            while (current.Parent != null)
            {
                height = height + 1;
                current = current.Parent;
            }

            return height;
        }
    }

}






