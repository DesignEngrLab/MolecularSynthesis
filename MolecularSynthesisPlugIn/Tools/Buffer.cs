using System;
using System.Collections.Generic;
using Priority_Queue;
using System.Diagnostics;
using System.IO;
using System.Linq;


namespace MolecularSynthesis.Tools
{
    public class JobBuffer
    {
        private const int MAX_SIMULATION = 150;
        private readonly SimplePriorityQueue<string, double> buffer;
        private readonly string _bufferDir;
        private Tuple<string, HashSet<string>>[] onSimulationTuples;
        private readonly Dictionary<string,int> epochLookUp;
        
        public JobBuffer(string runDir)
        {

            buffer = new SimplePriorityQueue<string, double>();
            _bufferDir = Path.Combine(runDir, "data");
            if (Directory.Exists(_bufferDir))
                Directory.Delete(_bufferDir, true);
            Directory.CreateDirectory(_bufferDir);
            onSimulationTuples = new Tuple<string, HashSet<string>>[]
            {
                Tuple.Create("short", new HashSet<string>()), 
                Tuple.Create("greaneylab", new HashSet<string>()),
            };

            epochLookUp = new Dictionary<string, int>();
        }

        public void Add(string linkerName, double priority, int epoch)
        {
            buffer.Enqueue(linkerName, priority);
            epochLookUp[linkerName] = epoch;
        }

        public bool Check_finised(Computation computation, StreamWriter sw, MessageClient clt)
        {
            foreach (var onSimulationInfo in onSimulationTuples)
            {
                var finished_linkers = new HashSet<string>();
                var queue = onSimulationInfo.Item1;
                var set = onSimulationInfo.Item2;
                foreach (var linkerName in set)
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
                        set.Remove(linkerName);
                        var property = computation.CalculateProperty(linkerName);
                        clt.SendMessage("[AddData]" + " " + linkerName);
                        Console.WriteLine("linker " + linkerName + " finished, with property " + property);
                        Console.WriteLine("Current "  + queue + " on simulation " + set.Count);
                        sw.WriteLine("Epoch " + epochLookUp[linkerName] + "," + linkerName + "," + property);
                        epochLookUp.Remove(linkerName);
                    }
                }
            }
            return Num_simulating() == 0;
        }

        public int Num_simulating()
        {
            return onSimulationTuples.Select(x => x.Item2.Count).ToArray().Sum();
        }

        public void Simulate()
        {
            var priority = buffer.GetPriority(buffer.First);
            var linkerName = buffer.Dequeue();
            var current_simulations = onSimulationTuples.Select(x => x.Item2.Count).ToArray();
            var target = Array.IndexOf(current_simulations, current_simulations.Min());

            var target_queue = onSimulationTuples[target].Item1;
            var target_set = onSimulationTuples[target].Item2;
            target_set.Add(linkerName);
            Submitlammps(linkerName, target_queue);

            Console.WriteLine("Job " + linkerName + " Submmitted with Priority " + priority);
            Console.WriteLine("Current "  + target_queue + " on simulation " + target_set.Count);
        }


        public bool CanFeedIn()
        {
            var current_simulations = onSimulationTuples.Select(x => x.Item2.Count);
            return buffer.Count > 0 && current_simulations.Min() < MAX_SIMULATION;
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
                //proc.WaitForExit();
            }

        }

    }


}