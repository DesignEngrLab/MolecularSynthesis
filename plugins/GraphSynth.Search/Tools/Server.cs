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
        private readonly int portUsed;
        

        public LearningServer(string learnDir, int port)
        {
            _learnDir = learnDir;
            portUsed = port;
            onlineSeverProcess = null;
        }


        public void StartOnlineServer()
        {
            if (onlineSeverProcess == null)
            {
                onlineSeverProcess = new Process();
                onlineSeverProcess.StartInfo.FileName = "/rhome/yangchen/.conda/envs/yangchenPython3/bin/python";
                onlineSeverProcess.StartInfo.Arguments = "server.py" + " " + portUsed.ToString();
                onlineSeverProcess.StartInfo.WorkingDirectory = Path.Combine(_learnDir);
                onlineSeverProcess.StartInfo.RedirectStandardError = true;
                onlineSeverProcess.StartInfo.UseShellExecute = false;
                onlineSeverProcess.StartInfo.RedirectStandardOutput = true;
                onlineSeverProcess.StartInfo.RedirectStandardInput = false;
                onlineSeverProcess.Start();
                //onlineSeverProcess.WaitForExit();
            }
            Console.WriteLine("Online server already started with Process ID: {0}", onlineSeverProcess.Id);
            //string error = onlineSeverProcess.StandardError.ReadToEnd();
            //Console.WriteLine(error);

            //Here the first output indicate using CPU or GPU
            Console.WriteLine(onlineSeverProcess.StandardOutput.ReadLine());
            //Here the first output indicate the time and host of the server, making sure the server has already started
            Console.WriteLine(onlineSeverProcess.StandardOutput.ReadLine());
        }

        public void MonitorOutput()
        {
            while (onlineSeverProcess != null)
            {
                var line = onlineSeverProcess.StandardOutput.ReadLine();
                if (line != "\n")
                    Console.WriteLine(line.Length);
            }
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