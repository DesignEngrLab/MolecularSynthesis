using System;
using System.Collections.Generic;
using System.IO;
using GraphSynth.Representation;
using OpenBabelFunctions;

namespace GraphSynth.Search.Algorithms
{
    public class BFS : AbstractAlgorithm
    {
        protected override string RunDirectoryName => "BFS";
        private Queue<BFSNode> QueueBFS;
        private HashSet<string> allNode;
        private HashSet<string> allFinalCand;
        private int nodeCnt;
        private int candCnt;
        private StreamWriter nodeInfoWriter;
        private StreamWriter candInfoWriter;

        public BFS(GlobalSettings settings) : base(settings)
        {
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

        public void search(candidate seed, int maxDepth)
        {
            BFSNode seedNode = new BFSNode(this, seed, 0);
            QueueBFS.Enqueue(seedNode);
            while (QueueBFS.Count != 0)
            {
                var popNode = QueueBFS.Dequeue();
                Console.WriteLine("{0},{1}", popNode.getSMILE(), popNode.getDepth());
//                Console.Write(".");
                if (popNode.getDepth() < maxDepth)
                {
                    var newNodes = popNode.getChildren();
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

                            var terminalLinker = node.getFinalCand();
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
        
        
    }
    

    public class BFSNode
    {
        private readonly BFS _parent;
        private candidate Cand;
        private int Depth;

        public BFSNode(BFS parent, candidate cand, int depth)
        {
            _parent = parent;
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

        public BFSNode[] getChildren()
        {
            var opts = _parent.GetAvailableOptions(Cand);
            var nodes = new BFSNode[opts.Count];
            for (var i = 0; i < nodes.Length; i++)
            {
                var child = _parent.CopyAndApplyOption(opts[i], Cand, true);
                nodes[i] = new BFSNode(_parent, child, Depth + 1);
            }
            return nodes;
        }

        public candidate[] getFinalCand()
        {
            var opts = _parent.GetCarboxylOptions(Cand);
            var finals = new candidate[opts.Count];
            for (var i = 0; i < finals.Length; i++)
            {
                finals[i] = _parent.CopyAndApplyOption(opts[i], Cand, true);
            }
            return finals;
        }
    }
}