using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;


namespace GraphSynth.Search.Tools
{
    public class LearningServer
    {
        private readonly string _runDir;
        private static readonly Dictionary<string,string> scriptLookup = new Dictionary<string, string>()
        {
            {"point", "calcPoint.py"}
        };
        private readonly string _learnDir;

        public LearningServer(string runDir, string learnDir)
        {
            _runDir = runDir;
            _learnDir = learnDir;
        }


        public void CalculateFeature(string task, string linkerId)
        {
            var _featureDir = Path.Combine(_runDir, task);
            if (Directory.Exists(_featureDir))
                Directory.Delete(_featureDir, true);
            Directory.CreateDirectory(_featureDir);
            
            var lmpData = Path.Combine(_runDir, "data", "linker" + linkerId + ".lmpdat");
            using (Process proc = new Process())
            {
                proc.StartInfo.FileName = "/rhome/yangchen/.conda/envs/yangchenPython3/bin/python";
                proc.StartInfo.Arguments = scriptLookup[task] + " " + lmpData + _featureDir;
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