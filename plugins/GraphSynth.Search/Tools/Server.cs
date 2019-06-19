using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Threading;


namespace GraphSynth.Search.Tools
{
    public class LearningServer
    {

        private readonly string _learnDir;
        private readonly string _runDir;
        private Process onlineSeverProcess;
        private readonly int portUsed;
        

        public LearningServer(string learnDir, int port, string runDir)
        {
            _learnDir = learnDir;
            _runDir = runDir;
            portUsed = port;
            onlineSeverProcess = null;
        }


        public void StartOnlineServer()
        {
            if (onlineSeverProcess == null)
            {
                onlineSeverProcess = new Process();
                onlineSeverProcess.StartInfo.FileName = "/rhome/yangchen/.conda/envs/yangchenPython3/bin/python";
                onlineSeverProcess.StartInfo.Arguments = "server.py" + " " + portUsed.ToString() + " " + _runDir;
                onlineSeverProcess.StartInfo.WorkingDirectory = Path.Combine(_learnDir);
                onlineSeverProcess.StartInfo.RedirectStandardError = false;
                onlineSeverProcess.StartInfo.UseShellExecute = false;
                onlineSeverProcess.StartInfo.RedirectStandardOutput = false;
                onlineSeverProcess.StartInfo.RedirectStandardInput = false;
                onlineSeverProcess.Start();
                //onlineSeverProcess.WaitForExit();
            }
            Thread.Sleep(15000);
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