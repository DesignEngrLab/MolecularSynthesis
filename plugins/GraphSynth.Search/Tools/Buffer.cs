using System;
using System.Collections.Generic;
using Priority_Queue;
using System.Diagnostics;
using System.Linq;
using System.Runtime.Remoting.Messaging;


namespace GraphSynth.Search.Tools
{
    public class JobBuffer
    {
        private const int MAX_SIMULATION = 10;
        private readonly SimplePriorityQueue<string, double> buffer;
        private readonly string _bufferDir;
        private HashSet<string> onSimulation;

        public JobBuffer(string dir)
        {
            buffer = new SimplePriorityQueue<string, double>();
            _bufferDir = dir;
        }

        public void Add(string linkerName, double priority)
        {
            buffer.Enqueue(linkerName, priority);
        }

        public bool Remove()
        {
            var priority = buffer.GetPriority(buffer.First);
            var linkerName = buffer.Dequeue();
            onSimulation.Add(linkerName);
            submitlammps(linkerName, "short");
            Console.WriteLine("Job " + linkerName + " Submmitted with Priority " + priority + ". Current on simulation " + onSimulation);
            if (linkerName == "finish")
                return true;
            return false;

        }

        public bool canFeedIn()
        {
            return buffer.Count > 0 && onSimulation.Count <= MAX_SIMULATION;
        }
        
        private void submitlammps(string linkerId, string queue) {
            using (Process proc = new Process()) {
                proc.StartInfo.FileName = "submit_lammps_linker_deform_remote";
                proc.StartInfo.Arguments = linkerId + " " + queue;
                proc.StartInfo.WorkingDirectory = _bufferDir;
                proc.StartInfo.RedirectStandardError = false;
                proc.StartInfo.UseShellExecute = false;
                proc.StartInfo.RedirectStandardOutput = true;
                proc.StartInfo.RedirectStandardInput = false;
                proc.Start();
                proc.WaitForExit();
            }

        }

    }


}