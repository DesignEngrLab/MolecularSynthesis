using System;
using System.Collections.Generic;
using Priority_Queue;
using System.IO;
using GraphSynth.Representation;
using GraphSynth.Search.Algorithms;
using OpenBabelFunctions;

namespace GraphSynth.Search.Tools
{
    public class JobBuffer
    {
        private const int MAX_SIMULATION = 20;
        private readonly SimplePriorityQueue<string, double> buffer;
        private int onSimulation;

        public JobBuffer()
        {
            buffer = new SimplePriorityQueue<string, double>();
        }

        public void Add(string linkerName, double priority)
        {
            buffer.Enqueue(linkerName, priority);
        }

        public string Remove()
        {
            var linkerName = buffer.Dequeue();
            onSimulation++;
            return linkerName;
        }

        public bool canFeedIn()
        {
            return buffer.Count > 0 && onSimulation < MAX_SIMULATION;
        }

    }


}