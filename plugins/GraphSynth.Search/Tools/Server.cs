using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;


namespace GraphSynth.Search.Tools
{
    public class LearningServer
    {
        private readonly string _runDir;
        private readonly string _learnDir;
        private readonly string _featureDir;
        private readonly string _task;
        private static readonly Dictionary<string,string> scriptLookup = new Dictionary<string, string>()
        {
            {"point", "calcPoint.py"}
        };

        public LearningServer(string runDir, string learnDir, string task)
        {
            _runDir = runDir;
            _learnDir = learnDir;
            _task = task;
            _featureDir = Path.Combine(_runDir, task);
            if (Directory.Exists(_featureDir))
                Directory.Delete(_featureDir, true);
            Directory.CreateDirectory(_featureDir);
        }


        public void CalculateFeature(string linkerId)
        {
            var lmpData = Path.Combine(_runDir, "data", "linker" + linkerId + ".lmpdat");
            using (Process proc = new Process())
            {
                proc.StartInfo.FileName = "/rhome/yangchen/.conda/envs/yangchenPython3/bin/python";
                proc.StartInfo.Arguments = scriptLookup[_task] + " " + lmpData + _featureDir;
                proc.StartInfo.WorkingDirectory = Path.Combine(_learnDir, "computation");
                proc.StartInfo.RedirectStandardError = true;
                proc.StartInfo.UseShellExecute = false;
                proc.StartInfo.RedirectStandardOutput = true;
                proc.StartInfo.RedirectStandardInput = false;
                proc.Start();
                proc.WaitForExit();
                string output = proc.StandardOutput.ReadToEnd();
                string error = proc.StandardError.ReadToEnd();
                Console.WriteLine(error);
                Console.WriteLine(output);
            }
        }

    }


}