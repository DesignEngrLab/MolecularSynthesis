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
using System.Diagnostics;
using System.Timers;
using System.Threading.Tasks;

namespace TestOpenBabel
{
    public class TestBenzRing : SearchProcess
    {
        public override string text => "TestBenz";
        public string filename = @"..\..\..\..\ForCiftest";
        public string extension = "xyz";
        //deault constructor
        public TestBenzRing(GlobalSettings settings) : base(settings)
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
            
            var stopwatch = new Stopwatch();

            TreeCandidate StartState = new TreeCandidate(seedCandidate);

            StartState.S = 0;
            StartState.n = 0;
            StartState.UCB = double.MaxValue;
            StartState.Children = new List<TreeCandidate>();

            var option0 = rulesets[0].recognize(StartState.graph);
            var option1 = rulesets[1].recognize(StartState.graph);
            var option2 = rulesets[2].recognize(StartState.graph);

            option0 = rulesets[0].recognize(StartState.graph);
            option0[4].apply(StartState.graph, null);
            //startstate.addtorecipe(option0[6]);

            option0 = rulesets[0].recognize(StartState.graph);
            option0[4].apply(StartState.graph, null);
            //startstate.addtorecipe(option0[5]);

            option0 = rulesets[0].recognize(StartState.graph);
            option0[2].apply(StartState.graph, null);
            //startstate.addtorecipe(option0[1]);

            option0 = rulesets[0].recognize(StartState.graph);
            option0[2].apply(StartState.graph, null);
            //startstate.addToRecipe(option0[1]);

            //option1 = rulesets[1].recognize(StartState.graph);
            //option1[4].apply(StartState.graph, null);
            ////StartState.addToRecipe(option1[3]);

            //option1 = rulesets[1].recognize(StartState.graph);
            //option1[9].apply(StartState.graph, null);
            ////StartState.addToRecipe(option1[16]);

            //option1 = rulesets[1].recognize(StartState.graph);
            //option1[10].apply(StartState.graph, null);
            ////StartState.addToRecipe(option1[16]);            

            option2 = rulesets[2].recognize(StartState.graph);
            option2[0].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[2]);


            var resultMol = OBFunctions.designgraphtomol(StartState.graph);
            resultMol = OBFunctions.InterStepMinimize(resultMol);
            OBFunctions.updatepositions(StartState.graph, resultMol);

            SearchIO.addAndShowGraphWindow(StartState.graph);

            var result = Evaluation.FindLengthAndRadius(StartState.graph);
            //[0] is Length, [1] is Radius, unit is 
            SearchIO.output("Length is: " + result[0]);
            SearchIO.output("Radius is: " + result[1]);

            //var FinalResultMol = OBFunctions.designgraphtomol(StartState.graph);

            //var conv = new OBConversion();
            //conv.SetInAndOutFormats("pdb", "mol");

            //int i = 12345678;
            //string name = ".mol";
            //name = Convert.ToString(i) + name;
            //conv.WriteFile(FinalResultMol, Path.Combine("C:\\Users\\zhang\\Desktop", name));

            //string name2 = ".xyz";
            //name2 = Convert.ToString(i) + name2;

            //using (Process proc = new Process())
            //{
            //    //"C:\Program Files\OpenBabel-3.1.1\obabel.exe"
            //    proc.StartInfo.FileName = "C:\\Program Files\\OpenBabel-3.1.1\\obabel.exe";
            //    proc.StartInfo.Arguments = name + " -O " + name2;
            //    proc.StartInfo.WorkingDirectory = "C:\\Users\\zhang\\Desktop";
            //    //C:\Users\zhang\Desktop
            //    //proc.StartInfo.RedirectStandardError = true;
            //    //proc.StartInfo.UseShellExecute = false;
            //    proc.StartInfo.RedirectStandardOutput = true;
            //    //proc.StartInfo.RedirectStandardInput = false;

            //    Console.Write("starting Converting...");
            //    proc.Start();

            //    //minimizeOutput = proc.StandardOutput.ReadToEnd();
            //    proc.WaitForExit();
            //}

            // writing TXT file
            //List<string> list = new List<string>();
            //list.Add("Hi");
            //list.Add("XXX");
            //System.IO.File.WriteAllLines(@"C:\Users\zhang\source\repos\MolecularSynthesis\output\MCTSProcess.txt", list);

            //resultMol = OBFunctions.InterStepMinimize(resultMol);
            //OBFunctions.updatepositions(seedGraph, resultMol);
            //var FinalResultMol = OBFunctions.designgraphtomol(seedGraph);

            //var conv = new OBConversion();
            //conv.SetInAndOutFormats("pdb", "mol");
            //conv.WriteFile(FinalResultMol, Path.Combine("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\output", "Test102.mol"));



            //SearchIO.addAndShowGraphWindow(seedGraph);

            //var result = new double[2];
            //result=Evaluation.FindLengthAndRadius(seedGraph);
            ////[0] is Length, [1] is Radius, unit is 
            //SearchIO.output("Length is: "+ result[0]);
            //SearchIO.output("Radius is: "+ result[1]);

            //string[] PathOfXyz = Directory.GetFiles(@"C:\Users\zhang\desktop", "*.xyz");
            //string CifFilename = Path.GetFileName(PathOfXyz[0]);
            //foreach (var namexxx in PathOfXyz)
            //{
            //    Console.WriteLine(Path.GetFileName(namexxx));
            //}
            //Console.WriteLine(Path.GetFileName(PathOfXyz[0]));
            //Console.WriteLine(Path.GetFileName(PathOfXyz[1]));

            //RunBat("C:\\Users\\zhang\\source\\repos\\PoreBlazer\\Windows\\HKUST1\\run.bat");

            //C:\Users\zhang\source\repos\PoreBlazer\Windows\HKUST1\results.txt
            //string[] EvaluationResult = System.IO.File.ReadAllLines(@"C:\Users\zhang\source\repos\PoreBlazer\Windows\HKUST1\results.txt");
            //Console.WriteLine(EvaluationResult[174]);

            stopwatch.Restart();
            var elapsed = stopwatch.Elapsed;
            Console.WriteLine("completed in {0}", elapsed);

            //var conv = new OBConversion();
            //conv.SetInAndOutFormats("pdb", extension);
            //conv.WriteFile(resultMol, filename + "." + extension);
            //File.AppendText()
        }

        private void RunBat(string batPath)
        {
            Process pro = new Process();

            FileInfo file = new FileInfo(batPath);
            pro.StartInfo.WorkingDirectory = file.Directory.FullName;
            pro.StartInfo.FileName = batPath;
            pro.StartInfo.CreateNoWindow = false;
            pro.StartInfo.UseShellExecute = true;
            pro.Start();
            pro.WaitForExit();
        }
    }




}

