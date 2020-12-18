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

namespace TestOpenBabel
{
    public class TestWholeProcess : SearchProcess
    {
        public override string text => "TestWholeProcess";
        //public string filename = @"..\..\..\..\ForCiftest";
        //public string extension = "xyz";


        //deault constructor
        public TestWholeProcess(GlobalSettings settings) : base(settings)
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

            List<string> MCTSProcess = new List<string>();

            TreeCandidate StartState = new TreeCandidate(seedCandidate);

            StartState.S = 0;
            StartState.n = 0;
            StartState.UCB = double.MaxValue;
            StartState.Children = new List<TreeCandidate>();

            var option0 = rulesets[0].recognize(StartState.graph);
            var option1 = rulesets[1].recognize(StartState.graph);
            var option2 = rulesets[2].recognize(StartState.graph);

            option0 = rulesets[0].recognize(StartState.graph);
            option0[6].apply(StartState.graph, null);
            StartState.addToRecipe(option0[6]);

            option0 = rulesets[0].recognize(StartState.graph);
            option0[5].apply(StartState.graph, null);
            StartState.addToRecipe(option0[5]);

            option0 = rulesets[0].recognize(StartState.graph);
            option0[1].apply(StartState.graph, null);
            StartState.addToRecipe(option0[1]);

            option1 = rulesets[1].recognize(StartState.graph);
            option1[3].apply(StartState.graph, null);
            StartState.addToRecipe(option1[3]);

            option1 = rulesets[1].recognize(StartState.graph);
            option1[16].apply(StartState.graph, null);
            StartState.addToRecipe(option1[16]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[1].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[1]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[2].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[2]);

            option2 = rulesets[2].recognize(StartState.graph);
            option2[0].apply(StartState.graph, null);
            StartState.addToRecipe(option2[0]);

            //option1 = rulesets[1].recognize(StartState.graph);
            //option1[35].apply(StartState.graph, null);
            //StartState.addToRecipe(option1[35]);

            //option1 = rulesets[1].recognize(StartState.graph);
            //option1[25].apply(StartState.graph, null);
            //StartState.addToRecipe(option1[25]);

            var resultMol = OBFunctions.designgraphtomol(StartState.graph);
            resultMol = OBFunctions.InterStepMinimize(resultMol);
            OBFunctions.updatepositions(StartState.graph, resultMol);

            var FinalResultMol = OBFunctions.designgraphtomol(StartState.graph);

            var conv = new OBConversion();
            conv.SetInAndOutFormats("pdb", "mol");

            int i = 12345678;
            string name = ".mol";
            name = Convert.ToString(i) + name;

            // "C:\\Users\\zhang\\source\\repos\\tobacco_3.0\\edges"
            conv.WriteFile(FinalResultMol, Path.Combine("C:\\Users\\zhang\\Desktop", name));

            string name2 = ".xyz";
            name2 = Convert.ToString(i) + name2;

            // Generate .xyz file after minimization
            using (Process proc = new Process())
            {
                //"C:\Program Files\OpenBabel-3.1.1\obabel.exe"
                proc.StartInfo.FileName = "C:\\Program Files\\OpenBabel-3.1.1\\obabel.exe";
                proc.StartInfo.Arguments = name + " -O " + name2;
                proc.StartInfo.WorkingDirectory = "C:\\Users\\zhang\\Desktop";
                //C:\\Users\\zhang\\Desktop
                proc.StartInfo.RedirectStandardOutput = true;

                proc.Start();

                proc.WaitForExit();
            }

            // Generate .cif file for tabacco input using julia

            using (Process proc = new Process())
            {
                // C:\Users\zhang\AppData\Local\Programs\Julia 1.5.3\bin\\julia.exe
                //C: \\Users\\zhang\\AppData\\Local\\Programs\\Julia 1.5.3\\bin\\julia.exe
                proc.StartInfo.FileName = "C:\\Users\\zhang\\AppData\\Local\\Programs\\Julia 1.5.3\\bin\\julia.exe";
                //proc.StartInfo.Arguments = name + " -O " + name2;

                // C:\Users\zhang\source\repos\MolecularSynthesis\CIFGeneration.jl
                // C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\CIFGeneration.jl
                proc.StartInfo.Arguments = "C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\CIFGeneration.jl";
                proc.StartInfo.WorkingDirectory = "C:\\Users\\zhang\\source\\repos\\tobacco_3.0\\edges";
                //C:\\Users\\zhang\\Desktop
                proc.StartInfo.RedirectStandardOutput = true;

                proc.Start();

                proc.WaitForExit();
            }

            // run tabacco

            using (Process proc = new Process())
            {
                //"C:\Users\zhang\.julia\conda\3\python.exe"

                //C:\Windows\System32\cmd.exe
                proc.StartInfo.FileName = "C:\\Windows\\System32\\cmd.exe";
                //proc.StartInfo.Arguments = name + " -O " + name2;
                proc.StartInfo.WorkingDirectory = "C:\\Users\\zhang\\source\\repos\\tobacco_3.0";
                //C:\\Users\\zhang\\Desktop
                proc.StartInfo.RedirectStandardInput = true;
                proc.StartInfo.RedirectStandardOutput = true;

                proc.Start();
                using (var sw = proc.StandardInput)
                {
                    if (sw.BaseStream.CanWrite)
                    {
                        // Vital to activate Anaconda

                        //C:\XYZ\anaconda\Scripts
                        sw.WriteLine("C:\\XYZ\\anaconda\\Scripts\\activate.bat");
                        // Activate your environment
                        sw.WriteLine("activate python27");
                        // Any other commands you want to run
                        sw.WriteLine("set KERAS_BACKEND=tensorflow");
                        // run your script. You can also pass in arguments
                        sw.WriteLine("python tobacco.py");
                    }
                }
                proc.WaitForExit();
            }

            // delete cif file in edges folder

            //string PreMOFName = "";
            //string SufMOFName = ".cif";
            //string MOFName = PreMOFName + SufMOFName;
            //System.IO.File.Delete(@"C:\Users\zhang\source\repos\tabacco_3.0\MOFName");

            // take output .cif file from tabacco to evaluation in poreblazer 4.0


            stopWatch.Stop();
            // Get the elapsed time as a TimeSpan value.
            TimeSpan ts = stopWatch.Elapsed;
            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
            ts.Hours, ts.Minutes, ts.Seconds,
            ts.Milliseconds / 10);

            MCTSProcess.Add(elapsedTime);
            System.IO.File.WriteAllLines(@"C:\Users\zhang\source\repos\MolecularSynthesis\output\WholeProcessRecord.txt", MCTSProcess);

        }
    }




}

