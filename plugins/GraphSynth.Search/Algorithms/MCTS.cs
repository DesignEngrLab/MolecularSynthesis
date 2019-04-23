/*************************************************************************
 *     This Mcts file & class is part of the GraphSynth.
 *     BaseClasses Project which is the foundation of the GraphSynth Ap-
 *     plication. GraphSynth.BaseClasses is protected and copyright under 
 *     the MIT License.
 *     Copyright (c) 2011 Matthew Ira Campbell, PhD.
 *
 *     Permission is hereby granted, free of charge, to any person obtain-
 *     ing a copy of this software and associated documentation files 
 *     (the "Software"), to deal in the Software without restriction, incl-
 *     uding without limitation the rights to use, copy, modify, merge, 
 *     publish, distribute, sublicense, and/or sell copies of the Software, 
 *     and to permit persons to whom the Software is furnished to do so, 
 *     subject to the following conditions:
 *     
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
 *     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGE-
 *     MENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
 *     FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
 *     CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION 
 *     WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *     
 *     Please find further details and contact information on GraphSynth
 *     at http://www.GraphSynth.com.
 *************************************************************************/

using System;
using GraphSynth.Search.Algorithms.Bandits;
using GraphSynth.Representation;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using GraphSynth.Search.Algorithms.Evaluations;
using OpenBabelFunctions;

namespace GraphSynth.Search.Algorithms {
    /// <summary>
    /// An overload for the RCA class that uses MCTS to choose the next action.
    /// </summary>
    public class Mcts : AbstractAlgorithm {
        private const int MaxDepth = 15; // how many levels to recurse into the constructed tree
//        private const int NumTrials = 50; // how many root node pulls per action choice
        private const int MaxWidth = 1;  // maximum number of pulls to perform at any given node for an action

        private int nodeCnt = 0;

        private readonly AbstractEvaluation _evaluation;  // how we evaluate leaf nodes

        private StreamWriter sw;

        /// <inheritdoc />
        public Mcts(GlobalSettings settings) : base(settings) {
            _evaluation = new RolloutEvaluation(1, 0, this);
        }
        
        public void buildTree(candidate cand) {
            sw = new StreamWriter("/home/manion/Documents/CleanMORFHPCC/output/MCTS/nodeInfo.txt");
            sw.WriteLine("NodeNum,Depth,Reward,SmileString");
            var rootNode = new BanditNode(this, cand, 0);
            
            var NumTrials = 1000;
            Console.WriteLine("Running {0} trails, Possible #options: {1}", NumTrials, rootNode.Options.Count);
            for (var i = 0; i < NumTrials; i++) {
                RunTrial(rootNode, MaxDepth);
            }
            Console.WriteLine("");
            Console.WriteLine("Avg reward and #tries for each option:");
            rootNode.Bandit.PrintOptInfo();
            sw.Close();
//            saveTree(rootNode);
            

        }
        
        /// <inheritdoc />
        public option ChooseOption(candidate cand) {
            var rootNode = new BanditNode(this, cand, 0);
//            var NumTrials = rootNode.Options.Count;
            var NumTrials = 50;
            Console.WriteLine("Running {0} trails, Possible #options: {1}", NumTrials, rootNode.Options.Count);
            for (var i = 0; i < NumTrials; i++) {
                Console.Write(".");
                RunTrial(rootNode, MaxDepth);
            }
            Console.WriteLine("");
            Console.WriteLine("Avg reward and #tries for each option:");
            rootNode.Bandit.PrintOptInfo();

            return rootNode.Options[rootNode.Bandit.GetBestArm()]; // return index of best arm we've found
        }

        private double RunTrial(BanditNode node, int depth) {
            if (depth == 0) // leaf node
                return _evaluation.Evaluate(node.Cand);
 

            if (!node.Bandit.HasOption()) {
                return _evaluation.Evaluate(node.Cand);
            }
            
            var optionIndex = node.Bandit.SelectPullArm();
            double totalReward;
            
//            Console.WriteLine("Querry opt idx: {0}, No.Children: {1}", optionIndex, node.Children.Length);
            // If we reach max child nodes, then select randomly among children according to how much we've visited
            if (node.Children[optionIndex].Count >= MaxWidth) {
//                Console.WriteLine("Should never be here1.");
                var successors = node.Children[optionIndex].Keys.ToList();
                var selectedOption = successors[node.Multinomial(optionIndex)];;
                node.Children[optionIndex][selectedOption].Visits += 1;
                var successorNode = node.Children[optionIndex][selectedOption].Node;
                totalReward = successorNode.TransitionReward + RunTrial(successorNode, depth - 1);
            } else {
                // generate a new successor node
                var successorState = CopyAndApplyOption(node.Options[optionIndex], node.Cand, true);
//                var immediateReward = Evaluate(successorState) - node.AbsoluteReward; // how much better than last node?
                var immediateReward = 0 - node.AbsoluteReward; // how much better than last node?
                
                // If the successor state is already in node.Children
                if (node.Children[optionIndex].ContainsKey(successorState)) {
//                    Console.WriteLine("Should never be here2.");
                    var successorNode = node.Children[optionIndex][successorState].Node;
                    node.Children[optionIndex][successorState].Visits += 1; // mark that we've sampled
                    totalReward = immediateReward + RunTrial(successorNode, depth - 1);
                } else {
                    var successorNode = new BanditNode(this, successorState, immediateReward);
                    node.Children[optionIndex][successorState] = new BanditNode.NodeCountTuple(successorNode);
                    totalReward = immediateReward + _evaluation.Evaluate(successorState);//this evalutation tells how this state is POTENTIALLY good
                    var fileDir = "_runDirectory" + "/intermediateLinkers/linker" + nodeCnt;
                    Directory.CreateDirectory(fileDir);
                    Settings.filer.Save("_runDirectory" + "/intermediateLinkers/linker" + nodeCnt + "/linker" + nodeCnt + ".xml", successorState);
                    Console.WriteLine("Node{0}: depth: {1}, reward: {2}, smi: {3}", nodeCnt, MaxDepth - depth + 1, totalReward,
                        OBFunctions.moltoSMILES(OBFunctions.designgraphtomol(successorState.graph)));
                    sw.WriteLine("{0},{1},{2},{3}", nodeCnt, MaxDepth - depth + 1, totalReward,
                        OBFunctions.moltoSMILES(OBFunctions.designgraphtomol(successorState.graph)));
                    nodeCnt++;
                }
            }
            
            node.Bandit.Update(optionIndex, totalReward);
            return totalReward;
        }


    }

    /// <summary>
    /// Stores information on a state, the reward for reaching the state, the options available, and the bandit used.
    /// </summary>
    public class BanditNode {
        private readonly Mcts _parent;

        public readonly candidate Cand;
        public readonly List<option> Options;
        public readonly AbstractBandit Bandit;
        public readonly double TransitionReward;
        public readonly double AbsoluteReward;  // absolute score on fitness function

        /// <summary>
        /// Each action is associated with a dictionary that stores successor bandits/states.
        /// The key for each successor is the state. 
        /// The value is a tuple [n,c], where n is a BanditNode and c is the number of times that n has been sampled.
        /// </summary>
        public Dictionary<candidate, NodeCountTuple>[] Children;

        public class NodeCountTuple {
            public readonly BanditNode Node;
            public int Visits = 1; // should be initialized on the first visit

            public NodeCountTuple(BanditNode node) {
                Node = node;
            }
        }
        

        public BanditNode(Mcts parent, candidate cand, double transitionReward) {
            _parent = parent;  // TODO transfer old tree
            Cand = cand;
            TransitionReward = transitionReward;
//            AbsoluteReward = _parent.Evaluate(cand);
            AbsoluteReward = 0;
            Options = AbstractAlgorithm.GetAvailableOptions(cand);
            Children = new Dictionary<candidate, NodeCountTuple>[Options.Count];
            for(var i=0; i<Children.Length; i++)
                Children[i] = new Dictionary<candidate, NodeCountTuple>();
//            Bandit = new EGreedyBandit(Options.Count, .5);  // TODO only track best and how to get there? credit
            Bandit = new UCTBandit(Options.Count);  
        }

        /// <summary>
        /// Samples the multinomial for the weighted number of visits to each child of the specified option index.
        /// </summary>
        public int Multinomial(int index) {
            var counts = Children[index].Keys.Select(k => Children[index][k].Visits).ToList();
            var countSum = counts.Sum();
            // List of counts proportional to number of times each child was sampled
            var averagedCounts = counts.Select(c => (double) c / countSum).ToList();

            var randVal = _parent.Rand.NextDouble();
            double current = 0;
            for (var i = 0; i < averagedCounts.Count; i++) {
                current += averagedCounts[i];
                if (randVal <= current)
                    return i;
            }
            return 0;  // won't ever get here
        }
        

    }



}