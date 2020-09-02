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

namespace MolecularSynthesis.Plugin
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

        public class TreeCandidate : candidate
        {
            public TreeCandidate Parent; //as a field
            public List<TreeCandidate> Children;
            public double s;
            public double N;
            public double UCB;

            public TreeCandidate(candidate seedCandidate)
            {
                this.activeRuleSetIndex = seedCandidate.activeRuleSetIndex;
                this.age = seedCandidate.age;
                this.designParameters = seedCandidate.designParameters;
                this.GenerationStatus = seedCandidate.GenerationStatus;
                this.graph = seedCandidate.graph;
                this.graphFileName = seedCandidate.graphFileName;
                this.performanceParams = seedCandidate.performanceParams;
            }

            public double CalculateUcb(TreeCandidate child)
            {
                return child.s / (double)child.N + 2 * Math.Sqrt(Math.Log(Parent.N) / child.N);
            }
            public TreeCandidate SelectPromisingNode()
            {
                TreeCandidate bestChild = null;
                double bestUcb = double.MinValue;

                foreach (TreeCandidate child in Children)
                {
                    if (child.N == 0)
                        return child;

                    double Ucb = CalculateUcb(child);
                    if (Ucb > bestUcb)
                    {
                        bestUcb = Ucb;
                        bestChild = child;
                    }
                }

                return bestChild;
            }

            public void BackPropogation(List<TreeCandidate> path)
            {
                foreach (var treeCandidate in path)
                {
                    treeCandidate.N++;
                    treeCandidate.s++;

                }
            }

            public double Rollout(TreeCandidate candidate, List<TreeCandidate> path)
            {
                double score;
                int RS0 = 0;
                //var childrenCandidate = RecognizeChooseApply.GenerateAllNeighbors(current, rulesets, false, false, true);

                foreach (var treecandidate in path)
                {
                    // need to recognize how many Rules from Ruleset0 exist
                    RS0 = RS0 + 1;

                }

                while (RS0 != 5)
                {
                    Random rnd = new Random();
                    //rnd.Next(0, 2); // generate 0 or 1
                    if (rnd.Next(0, 2) == 0)
                    {
                        RS0 = RS0 + 1;
                    }
                    //candidate=RecognizeChooseApply.GenerateAllNeighbors(current, rulesets, false, false, true)
                }

                score = Evaluation.distance(candidate, desiredLenghtAndRadius);
                return score;
            }

            public override candidate copy()
            {
                return new TreeCandidate(base.copy());
            }
        }
        protected override void Run()
        {
            //var candidates = new SimplePriorityQueue<candidate, double>();

            // generate a random number 0 or 1 to decide the next rule is from RS0 or RS1
            Random rnd = new Random();
            //rnd.Next(0, 2); // generate 0 or 1

            // use 10000 is that DS use 3000-70000 iteration for 9*9 go play , so guess 10000 is enough
            int iteration = 10000;
            //TreeCandidate node1 = new TreeCandidate() { S = 0, n=0, UCB=0 };


            // 1. check if this is the leaf node, if no go to step 2 until it is a leaf node,if yes go to step 3
            // 2. find the children who has the best UCB value
            // 3. do random simulation
            // 4. update S,n,UCB value for the whole tree
            TreeCandidate current = (TreeCandidate)seedCandidate;
            current.s = 0;
            current.N = 0;
            current.UCB = int.MaxValue;
            current.Children = null;

            for (int i = 0; i < iteration; i++)
            {
                current.SelectPromisingNode();
                current.Rollout(current, current.Children);
                current.BackPropogation(current.Children);

            }

        }
    }
}