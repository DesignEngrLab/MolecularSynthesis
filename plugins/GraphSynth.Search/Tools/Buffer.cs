using System.Collections.Generic;
using Priority_Queue;
using System.IO;
using GraphSynth.Representation;
using GraphSynth.Search.Algorithms;
using OpenBabelFunctions;

namespace GraphSynth.Search.Tools
{
    public class MolBuffer
    {
        private readonly string _fileDir;
        private readonly SimplePriorityQueue<string, double> buffer;
        private readonly Dictionary<string, candidate> candLookUpDict;

        public MolBuffer(string fileDir)
        {
            _fileDir = fileDir;
            buffer = new SimplePriorityQueue<string, double>();
            candLookUpDict = new Dictionary<string, candidate>();
        }

        public void Add(string linkerName, candidate cand, double priority)
        {
            buffer.Enqueue(linkerName, priority);
            candLookUpDict.Add(linkerName, cand);
        }

        public void Remove()
        {
            var linkerName = buffer.Dequeue();
            var cand = candLookUpDict[linkerName];
            WriteFile(linkerName, cand);
            candLookUpDict.Remove(linkerName);
        }

        public int Size()
        {
            return candLookUpDict.Count;
        }
        
        private void WriteFile(string linkerName, candidate cand)
        {
            var coeff = Path.Combine(_fileDir, "linker-" + linkerName + ".coeff");
            var lmpdat = Path.Combine(_fileDir, "linker-" + linkerName + ".lmpdat");
            AbstractAlgorithm.Converter.moltoUFF(OBFunctions.designgraphtomol(cand.graph), coeff, lmpdat, false, 100);
        }
    }


}