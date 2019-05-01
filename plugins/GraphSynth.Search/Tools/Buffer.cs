using System;
using System.Collections.Generic;
using Priority_Queue;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.Remoting.Messaging;


namespace GraphSynth.Search.Tools
{
    public class JobBuffer
    {
        private const int MAX_SIMULATION = 10;
        private readonly SimplePriorityQueue<string, double> buffer;
        private readonly string _bufferDir;
        private readonly HashSet<string> onSimulation;

        public JobBuffer(string dir)
        {
            buffer = new SimplePriorityQueue<string, double>();
            _bufferDir = dir;
            onSimulation = new HashSet<string>();
        }

        public void Add(string linkerName, double priority)
        {
            buffer.Enqueue(linkerName, priority);
        }

        public void Check_finised()
        {
            foreach (var linkerName in onSimulation)
            {
                var simulationDir = Path.Combine(_bufferDir, "linker" + linkerName + "_deformation");
                if (File.Exists(Path.Combine(simulationDir, "DONE")))
                {
                    Console.WriteLine("linker" + linkerName + "finished");
                    onSimulation.Remove(linkerName);
                }
            }
        }

        public bool Remove()
        {
            var priority = buffer.GetPriority(buffer.First);
            var linkerName = buffer.Dequeue();
            onSimulation.Add(linkerName);
            Submitlammps(linkerName, "short");
            Console.WriteLine("Job " + linkerName + " Submmitted with Priority " + priority + ". Current on simulation " + onSimulation.Count);
            if (linkerName == "finish")
                return true;
            return false;

        }

        public bool CanFeedIn()
        {
            return buffer.Count > 0 && onSimulation.Count < MAX_SIMULATION;
        }
        
        private void Submitlammps(string linkerId, string queue) {
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