using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;


namespace GraphSynth.Search.Tools
{
    public class LearningServer
    {

        private readonly string _learnDir;
        private Process onlineSeverProcess;
        

        public LearningServer(string learnDir)
        {
            _learnDir = learnDir;
        }


        public void StartOnlineServer()
        {
            if (onlineSeverProcess == null)
            {
                onlineSeverProcess = new Process();
                onlineSeverProcess.StartInfo.FileName = "/rhome/yangchen/.conda/envs/yangchenPython3/bin/python";
                onlineSeverProcess.StartInfo.Arguments = "server.py";
                onlineSeverProcess.StartInfo.WorkingDirectory = Path.Combine(_learnDir);
                onlineSeverProcess.StartInfo.RedirectStandardError = true;
                onlineSeverProcess.StartInfo.UseShellExecute = false;
                onlineSeverProcess.StartInfo.RedirectStandardOutput = true;
                onlineSeverProcess.StartInfo.RedirectStandardInput = false;
                onlineSeverProcess.Start();
                //proc.WaitForExit();
                //string error = proc.StandardError.ReadToEnd();
                //Console.WriteLine(error);
                //string output = proc.StandardOutput.ReadToEnd();
                //Console.WriteLine(output);
            }
            Console.WriteLine("Online server already started with Process ID: {0}", onlineSeverProcess.Id);
        }

        public void ShutDownOnlineServer()
        {
            if (onlineSeverProcess != null)
                onlineSeverProcess.Kill();
            Console.WriteLine("Online server with Process ID: {0} has already shutted down", onlineSeverProcess.Id);
            onlineSeverProcess = null;
        }
    }


}