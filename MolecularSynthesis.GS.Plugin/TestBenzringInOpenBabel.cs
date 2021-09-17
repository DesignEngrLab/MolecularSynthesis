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


        static double[] desiredLenghtAndRadius = new double[] { 565, 140 };
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

            //3 4 5 6 (658, 94)


            // 3 6 3 3 4
            // 2 5 2 2 3
            // 638.247 93.74


            //3 4 3 3 6 
            //2 3 2 2 5 

            //3 1 1 3 1
            //2 0 0 2 0

            //3 5 3 4 3(644 93)
            //2 4 2 3 2

            // 3 4 2 2 1 (601,81)
            // 2 3 1 1 0

            // 6 3 6 3 2(602,68)
            // 5 2 5 2 1

            // 3 3 5 3 6

            // 2 2 4 2 5

            //RS0 7 3 1 5 2(565,140)  RS1 option 4,9,10 (recegonize before apply everytime!!!)
            //    6 2 0 4 1 
            // RS0 2 5 2 2 7 --- RS1 4  (583.79, 111.58)
            //     1 4 1 1 6         12

            // RS0 7 7 5 5 ; RS1 8(17.33)
            //     6 6 4 4 


            // Rs0 3 7 4 6 ; RS1 7(option number 24)   (5.5 or 7)
            //     2 6 3 5       24    
            // 1 7 3 6 1-- 8 3 (D: 6.97)  (572,117)
            // 0 6 2 5 0   7 2


            // 6 5 7 1 ; 1_7 1_3
            option0 = rulesets[0].recognize(StartState.graph);
            option0[5].apply(StartState.graph, null);
            StartState.addToRecipe(option0[5]);

            option0 = rulesets[0].recognize(StartState.graph);
            option0[4].apply(StartState.graph, null);
            StartState.addToRecipe(option0[4]);

            option0 = rulesets[0].recognize(StartState.graph);
            option0[6].apply(StartState.graph, null);
            StartState.addToRecipe(option0[6]);

            option0 = rulesets[0].recognize(StartState.graph);
            option0[0].apply(StartState.graph, null);
            StartState.addToRecipe(option0[0]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[5].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[5]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[2].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[2]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[2].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[2]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[4].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[4]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[2].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[2]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[5].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[5]);

            Console.WriteLine("\n");
            Console.WriteLine("--------------------");

            var candidate = (TreeCandidate)StartState.copy();

            //foreach (var item in candidate.recipe)
            //{
            //    // item.ruleSetIndex
            //    // item.ruleNumber
            //    Console.WriteLine("RulesetNumber: " + item.ruleSetIndex.ToString() +"    "+ "RuleNumber: "+ item.ruleNumber.ToString());



            //}
            //Console.WriteLine("--------------------");

            //option2 = rulesets[2].recognize(candidate.graph);
            //option2[0].apply(candidate.graph, null);
            //candidate.addToRecipe(option0[2]);

            option1 = rulesets[1].recognize(StartState.graph);
            option1[25].apply(StartState.graph, null);
            StartState.addToRecipe(option1[25]);

            //option1 = rulesets[1].recognize(StartState.graph);
            option1[7].apply(StartState.graph, null);
            StartState.addToRecipe(option1[7]);


            //var resultMol = OBFunctions.designgraphtomol(candidate.graph);
            //resultMol = OBFunctions.InterStepMinimize(resultMol);
            //OBFunctions.updatepositions(candidate.graph, resultMol);

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


            //----------------------------------------------------------------------
            option2 = rulesets[2].recognize(StartState.graph);
            option2[0].apply(StartState.graph, null);
            StartState.addToRecipe(option0[2]);


            var resultMol = OBFunctions.designgraphtomol(StartState.graph);
            resultMol = OBFunctions.InterStepMinimize(resultMol);
            OBFunctions.updatepositions(StartState.graph, resultMol);

            SearchIO.addAndShowGraphWindow(StartState.graph);

            var result = Evaluation.FindLengthAndRadius(StartState.graph);
            var score = Evaluation.distance(StartState, desiredLenghtAndRadius);
            //[0] is Length, [1] is Radius, unit is 
            SearchIO.output("Length is: " + result[0]);
            SearchIO.output("Radius is: " + result[1]);
            SearchIO.output("distance is: " + score);

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

