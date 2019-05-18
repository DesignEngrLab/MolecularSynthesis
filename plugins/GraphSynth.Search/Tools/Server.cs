using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.Remoting.Messaging;


namespace GraphSynth.Search.Tools
{
    public class LearningServer
    {
        private readonly string _runDir;
        private readonly string _featureDir;
        private readonly string _learnDir;

        public LearningServer(string runDir, string learnDir)
        {
            _runDir = runDir;
            _learnDir = learnDir;
        }


        public void CalculateFeature(string script, string linkerId)
        {
            var lmpData = Path.Combine(_runDir, "data", "linker" + linkerId + ".lmpdat");
            using (Process proc = new Process())
            {
                proc.StartInfo.FileName = "/rhome/yangchen/.conda/envs/yangchenPython3/bin/python";
                proc.StartInfo.Arguments = script + " " + lmpData;
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