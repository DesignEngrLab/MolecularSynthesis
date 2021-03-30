using System.Collections.Generic;
using System.Threading;
using GraphSynth;
using GraphSynth.Representation;
using GraphSynth.Search;
using Priority_Queue;
using System;
using System.Collections;
using System.Security.Cryptography.X509Certificates;
using System.IO;
using MolecularSynthesis.GS.Plugin;
using System.Linq;
using OpenBabel;
using OpenBabelFunctions;
using MolecularSynthesis.GS.Plugin;
using System.Diagnostics;
using System.Timers;
using System.Diagnostics;
using System.Threading.Tasks;
using System.Text;
using System.Reflection;

namespace TestOpenBabel
{
    public class Test_zeo : SearchProcess
    {
        public override string text => "Testzeo";
        //public string filename = @"..\..\..\..\ForCiftest";
        //public string extension = "xyz";


        //deault constructor
        public Test_zeo(GlobalSettings settings) : base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 2;
            AutoPlay = true;
        }
        protected override void Run()
        {
            //var resultMol = OBFunctions.designgraphtomol(seedGraph);
            // after convert from designgraph to .mol file
            // save the result as .mol file
            // call minimize.exe to do the energy minimization               

            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();


            var dirString = this.rulesDirectory;
            var dir = new DirectoryInfo(dirString);
            dir = dir.Parent;
            var zeoDir = new DirectoryInfo(Path.Combine(dir.FullName, "zeo++-0.3")).FullName;

            var filepath= Path.Combine(zeoDir, "network");
            var arguments= "-res " + Path.Combine(zeoDir, "IRMOF-1.cssr");
            Console.WriteLine(filepath);
            Console.WriteLine(arguments);

            using (Process proc = new Process())
            {
                // "C:\Users\zhang\AppData\Local\Programs\Julia 1.5.3\bin\\julia.exe"
                //C: \\Users\\zhang\\AppData\\Local\\Programs\\Julia 1.5.3\\bin\\julia.exe

                //C: \Users\zhang\source\repos\zeo++-0.3\network
                proc.StartInfo.FileName =Path.Combine(zeoDir, "network");
                //proc.StartInfo.Arguments = name + " -O " + name2;

                // "C:\Users\zhang\source\repos\MolecularSynthesis\CIFGeneration.jl"
                // C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\CIFGeneration.jl
                proc.StartInfo.Arguments = "-res " + Path.Combine(zeoDir, "IRMOF-1.cssr");
                //C: \Users\zhang\source\repos\MolecularSynthesis\output
                proc.StartInfo.WorkingDirectory = this.outputDirectory;
                Console.WriteLine(this.outputDirectory);
                //C:\\Users\\zhang\\Desktop
                proc.StartInfo.RedirectStandardOutput = true;
                proc.StartInfo.UseShellExecute=false;
                proc.Start();

                proc.WaitForExit();
            }



            stopWatch.Stop();
            // Get the elapsed time as a TimeSpan value.
            TimeSpan ts = stopWatch.Elapsed;
            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
            ts.Hours, ts.Minutes, ts.Seconds,
            ts.Milliseconds / 10);

            Console.WriteLine("Total Running Time" + elapsedTime);




        }
    }
}

