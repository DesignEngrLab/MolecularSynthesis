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
        private const int MAX_SIMULATION = 150;
        private readonly SimplePriorityQueue<string, double> buffer;
        private readonly string _bufferDir;
        private readonly HashSet<string> onSimulation;


        public JobBuffer(string runDir)
        {
            buffer = new SimplePriorityQueue<string, double>();
            _bufferDir = Path.Combine(runDir, "data");
            if (Directory.Exists(_bufferDir))
                Directory.Delete(_bufferDir, true);
            Directory.CreateDirectory(_bufferDir);
            onSimulation = new HashSet<string>();
        }

        public void Add(string linkerName, double priority)
        {
            buffer.Enqueue(linkerName, priority);
        }

        public bool Check_finised(LearningServer server)
        {
            var finished_linkers = new HashSet<string>();
            foreach (var linkerName in onSimulation)
            {
                var simulationDir = Path.Combine(_bufferDir, "linker" + linkerName + "_deformation");
                if (File.Exists(Path.Combine(simulationDir, "DONE")))
                {
                    finished_linkers.Add(linkerName);
                }
            }
            if (finished_linkers.Count > 0)
            {
                foreach (var linkerName in finished_linkers)
                {
                    onSimulation.Remove(linkerName);
                    var property = server.CalculateProperty(linkerName);
                    Console.WriteLine("linker " + linkerName + " finished, with property " + property);
                    Console.WriteLine("Current on simulation " + onSimulation.Count);
                }
            }
            return onSimulation.Count == 0;
        }

        public bool Simulate()
        {
            var priority = buffer.GetPriority(buffer.First);
            var linkerName = buffer.Dequeue();
            if (linkerName == "finish")
                return true;
            onSimulation.Add(linkerName);
            string[] epochAndName = linkerName.Split(",");
            Submitlammps(epochAndName[1], "short");
            //Console.WriteLine("Job " + linkerName + " Submmitted with Priority " + priority);
            //Console.WriteLine("Current on simulation " + onSimulation.Count);
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