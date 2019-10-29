﻿using System;
using System.Collections.Generic;
using System.IO;
using GraphSynth.Representation;
using GraphSynth.Search.Algorithms;
using OpenBabelFunctions;


namespace GraphSynth.Search
{
    public class Bfs: SearchProcess
    {
        private readonly string _runDirectory;
        private readonly string _inputFilePath;
        
        private candidate Seed;
        private const int MAX_DEPTH = 2;
        private Queue<BFSNode> QueueBFS;
        private HashSet<string> allNode;
        private HashSet<string> allFinalCand;
        private int nodeCnt;
        private int candCnt;
        private StreamWriter nodeInfoWriter;
        private StreamWriter candInfoWriter;
        
        
        /// <inheritdoc />
        /// <summary>
        /// Initializes SearchProcess properties.
        /// </summary>
        public Bfs(GlobalSettings settings): base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 1;
            AutoPlay = true;
            
            _runDirectory = Path.Combine(settings.OutputDirAbs, "BFS");
            if (!Directory.Exists(_runDirectory))
                Directory.CreateDirectory(_runDirectory);
            Seed = new candidate(OBFunctions.tagconvexhullpoints(settings.seed), settings.numOfRuleSets);
            QueueBFS = new Queue<BFSNode>();
            allNode = new HashSet<string>();
            allFinalCand = new HashSet<string>();
            nodeCnt = 0;
            candCnt = 0;
            nodeInfoWriter = new StreamWriter(_runDirectory + "/nodeInfo.txt");
            nodeInfoWriter.WriteLine("SMILE,Depth");
            candInfoWriter = new StreamWriter(_runDirectory + "/candInfo.txt");
            candInfoWriter.WriteLine("SMILE,Depth");
        }

        protected override void Run()
        {
            var agent = new Algorithms.Deterministic(settings);
            BFSNode seedNode = new BFSNode(Seed, 0);
            QueueBFS.Enqueue(seedNode);
            while (QueueBFS.Count != 0)
            {
                var popNode = QueueBFS.Dequeue();
                Console.WriteLine("{0},{1}", popNode.getSMILE(), popNode.getDepth());
//                Console.Write(".");
                if (popNode.getDepth() < MAX_DEPTH)
                {
                    var newNodes = popNode.getChildren(agent);
                    foreach (var node in newNodes)
                    {
                        var nodeSMILE = node.getSMILE();
                        if (!allNode.Contains(nodeSMILE))
                        {
                            nodeInfoWriter.WriteLine("{0},{1}", nodeSMILE, node.getDepth());
//                            var nodeDir = _runDirectory + "/interNode/depth" + node.getDepth() + "/node" + nodeCnt;
//                            Directory.CreateDirectory(nodeDir);
//                            Settings.filer.Save(nodeDir + "/node" + nodeCnt + ".xml", node.getCand());
                            allNode.Add(nodeSMILE);
                            QueueBFS.Enqueue(node);
                            nodeCnt++;

                            var terminalLinker = node.getFinalCand(agent);
                            foreach (var linker in terminalLinker)
                            {
                                var mol = OBFunctions.designgraphtomol(linker.graph);
                                var candSMILE = OBFunctions.moltoSMILES(mol);
                                if (!allFinalCand.Contains(candSMILE))
                                {
//                                    var linkerDir = _runDirectory + "/finalNode/depth" + node.getDepth() + "/linker" + candCnt;
//                                    Directory.CreateDirectory(linkerDir);
//                                    
//                                    var coeff = Path.Combine(linkerDir, "linker" + candCnt + ".coeff");
//                                    var lmpdat = Path.Combine(linkerDir, "linker" + candCnt + ".lmpdat");
                                    allFinalCand.Add(candSMILE);
//                                    // Set up UFF and run lammps
//                                    Converter.moltoUFF(mol, coeff, lmpdat, false, 100);
//                                    Settings.filer.Save(linkerDir + "/linker" + candCnt + ".xml", linker.graph);
                                    candInfoWriter.WriteLine("{0},{1}", candSMILE, node.getDepth());
                                    candCnt++;
                                }
                            }

                        }

                    }
                }
            }
            nodeInfoWriter.Close();
            candInfoWriter.Close();
            
            
        }

        public override string text => "BFS Search Runner";
    }
    
    
    public class BFSNode
    {
        private candidate Cand;
        private int Depth;

        public BFSNode(candidate cand, int depth)
        {
            Cand = cand;
            Depth = depth;
        }

        public int getDepth()
        {
            return Depth;
        }

        public candidate getCand()
        {
            return Cand;
        }

        public string getSMILE()
        {
            return OBFunctions.moltoSMILES(OBFunctions.designgraphtomol(Cand.graph));
        }

        public BFSNode[] getChildren(Deterministic agent)
        {
            var opts = AbstractAlgorithm.GetNoneTerminalOptions(Cand);
            var nodes = new BFSNode[opts.Count];
            for (var i = 0; i < nodes.Length; i++)
            {
                var child = agent.CopyAndApplyOption(opts[i], Cand, true);
                nodes[i] = new BFSNode(child, Depth + 1);
            }
            return nodes;
        }

        public candidate[] getFinalCand(Deterministic agent)
        {
            var opts = AbstractAlgorithm.GetTerminalOptions(Cand);
            var finals = new candidate[opts.Count];
            for (var i = 0; i < finals.Length; i++)
            {
                finals[i] = agent.CopyAndApplyOption(opts[i], Cand, true);
            }
            return finals;
        }
    }
}