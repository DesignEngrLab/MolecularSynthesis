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
    public class Test_zeo_as_evaluation_para : SearchProcess
    {
        public override string text => "Test_zeo_as_evaluation_para";
        //public string filename = @"..\..\..\..\ForCiftest";
        //public string extension = "xyz";
        static Random rand = new Random();
        static TreeCandidate noneparallel = new TreeCandidate(new candidate());

        //deault constructor
        public Test_zeo_as_evaluation_para(GlobalSettings settings) : base(settings)
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
            List<string> Results = new List<string>();

            Parallel.For(0, 10, count =>
            {

                TreeCandidate StartState = new TreeCandidate(seedCandidate);

                StartState.S = 0;
                StartState.n = 0;
                StartState.UCB = double.MaxValue;
                StartState.Children = new List<TreeCandidate>();

                var option0 = rulesets[0].recognize(StartState.graph);
                var option1 = rulesets[1].recognize(StartState.graph);
                var option2 = rulesets[2].recognize(StartState.graph);

                var candidate = (TreeCandidate)StartState.copy();

                int ThreadNumber = System.Threading.Thread.CurrentThread.ManagedThreadId;
                candidate = (TreeCandidate)StartState.copy();
                option0 = rulesets[0].recognize(candidate.graph);
                option1 = rulesets[1].recognize(candidate.graph);
                option2 = rulesets[2].recognize(candidate.graph);

                //List<string> Results = new List<string>();

                for (int j = 0; j < 4; j++)
                {


                    //rnd.Next(0, 2); 0 or 1
                    var RuleSetNumber = rand.Next(0, 1);
                    var TotalOption = rulesets[RuleSetNumber].recognize(candidate.graph).Count;
                    var OptionNumber = rand.Next(0, TotalOption);
                    rulesets[RuleSetNumber].recognize(candidate.graph)[OptionNumber].apply(candidate.graph, null);
                    candidate.addToRecipe(rulesets[RuleSetNumber].recognize(candidate.graph)[OptionNumber]);

                }

                option2 = rulesets[2].recognize(candidate.graph);
                option2[0].apply(candidate.graph, null);
                //StartState.addToRecipe(option2[0]);

                var resultMol = OBFunctions.designgraphtomol(candidate.graph);
                resultMol = justMinimize(resultMol);
                OBFunctions.updatepositions(candidate.graph, resultMol);

                var FinalResultMol = OBFunctions.designgraphtomol(candidate.graph);

                var conv = new OBConversion();
                conv.SetInAndOutFormats("pdb", "mol");

                // 1. generate .mol file, move the zeo++ folder
                int i = ThreadNumber;
                string name = ".mol";
                name = Convert.ToString(i) + name;
                conv.WriteFile(FinalResultMol, Path.Combine("/nfs/hpc/share/zhangho2/MolecularSynthesis/output", name));

                string position1 = "/nfs/hpc/share/zhangho2/MolecularSynthesis/output/" + name;
                string position2 = "/nfs/hpc/share/zhangho2/zeo++-0.3/" + name;
                System.IO.File.Move(position1, position2, true);
                Console.WriteLine("\n");
                //2. .mol to.xyz

                string name2 = ".xyz";
                name2 = Convert.ToString(i) + name2;
                string position3 = "/nfs/hpc/share/zhangho2/zeo++-0.3/" + name2;

                using (Process proc = new Process())
                {
                    //"C:\Program Files\OpenBabel-3.1.1\obabel.exe"
                    proc.StartInfo.FileName = "/usr/local/apps/openbabel/3.1.1/bin/obabel";
                    proc.StartInfo.Arguments = position2 + " -O " + position3;
                    proc.StartInfo.WorkingDirectory = "/nfs/hpc/share/zhangho2/zeo++-0.3";
                    //C:\Users\zhang\Desktop
                    //proc.StartInfo.RedirectStandardError = true;
                    //proc.StartInfo.UseShellExecute = false;
                    proc.StartInfo.RedirectStandardOutput = true;
                    //proc.StartInfo.RedirectStandardInput = false;

                    Console.Write("starting Converting...");
                    proc.Start();

                    //minimizeOutput = proc.StandardOutput.ReadToEnd();
                    proc.WaitForExit();
                }

                Console.WriteLine("\n");

                // 3. get rid of two carboxylate
                string name3 = Convert.ToString(i) + "_XXX" + ".xyz";
                string position4 = "/nfs/hpc/share/zhangho2/zeo++-0.3/" + name3;

                using (Process proc = new Process())
                {
                    proc.StartInfo.FileName = "/nfs/hpc/share/zhangho2/zeo++-0.3/molecule_to_abstract";
                    proc.StartInfo.Arguments = position3 + " 0 " + position4;
                    proc.StartInfo.WorkingDirectory = this.outputDirectory;
                    proc.StartInfo.RedirectStandardOutput = true;
                    proc.Start();
                    Console.Write("starting removing...");
                    proc.WaitForExit();

                }
                Console.WriteLine("\n");
                // 4. read the new xyz file and copy

                string[] lines = System.IO.File.ReadAllLines(@position4);
                System.IO.File.WriteAllLines(@"/nfs/hpc/share/zhangho2/zeo++-0.3/ForEvaluation.xyz", lines);
                Console.WriteLine("Finish writing new .xyz file");
                Console.WriteLine("\n");

                // 5.  build MOF


                // ./framework_builder nets/pcu.cgd 1 output 6c_Zn_1_Ch.xyz ForEvaluation.xyz
                using (Process proc = new Process())
                {
                    proc.StartInfo.FileName = "/nfs/hpc/share/zhangho2/zeo++-0.3/framework_builder";
                    proc.StartInfo.Arguments = "nets/pcu.cgd 1 output 6c_Zn_1_Ch.xyz ForEvaluation.xyz";
                    proc.StartInfo.WorkingDirectory = "/nfs/hpc/share/zhangho2/zeo++-0.3";
                    proc.StartInfo.RedirectStandardOutput = true;
                    proc.StartInfo.UseShellExecute = false;
                    proc.StartInfo.RedirectStandardError = true;
                    proc.Start();
                    Console.Write("starting building...");
                    proc.WaitForExit();

                }

                //  5.1 need to change the output file name for multithread

                

                    string finalVar = ThreadNumber.ToString() + ".cssr";
                lock (noneparallel)
                    File.Delete("/nfs/hpc/share/zhangho2/zeo++-0.3/" + finalVar);
                    System.IO.File.Move("/nfs/hpc/share/zhangho2/zeo++-0.3/output_framework.cssr", "/nfs/hpc/share/zhangho2/zeo++-0.3/" + finalVar);




                // 6. find accessible volume 
                using (Process proc = new Process())
                {

                    //C: \Users\zhang\source\repos\zeo++-0.3\network
                    // /nfs/hpc/share/zhangho2/MolecularSynthesis/zeo++-0.3
                    proc.StartInfo.FileName = "/nfs/hpc/share/zhangho2/zeo++-0.3/network";
                    //proc.StartInfo.Arguments = name + " -O " + name2;

                    proc.StartInfo.Arguments = " -vol 1.2 1.2 50000 " + finalVar;
                    //C: \Users\zhang\source\repos\MolecularSynthesis\output
                    proc.StartInfo.WorkingDirectory = "/nfs/hpc/share/zhangho2/zeo++-0.3";
                    //C:\\Users\\zhang\\Desktop
                    proc.StartInfo.RedirectStandardOutput = true;
                    proc.StartInfo.UseShellExecute = false;
                    proc.Start();

                    proc.WaitForExit();
                }

                // 7. read data from relative file


                string contents = File.ReadAllText("/nfs/hpc/share/zhangho2/zeo++-0.3/" + ThreadNumber.ToString() + ".vol");
                string[] words = contents.Split(' ');
                Console.WriteLine(contents);

                Console.WriteLine("------------------");
                foreach (var word in words)
                {
                    Console.WriteLine(word);
                }

                Console.WriteLine("Accessible volume:" + words[15]);
                Console.WriteLine("Accessible Volume Fraction:  " + words[13]);
                Results.Add(words[15]);
                Results.Add(words[13]);
                ;            //Console.WriteLine("Poresize: ", words[5]);
                             //PoreSizeValue = Convert.ToDouble(words[5]);
                             //Console.WriteLine("Poresizevalue: ", words[5]);

                //Results.Add("Accessible volume: " + words[15] + "---" + ThreadNumber.ToString());
                //Results.Add("Accessible Volume Fraction: " + words[13] + "---" + ThreadNumber.ToString());


                File.Delete("/nfs/hpc/share/zhangho2/zeo++-0.3/" + finalVar + ".cssr");

            });





            System.IO.File.WriteAllLines(@"/nfs/hpc/share/zhangho2/zeo++-0.3/AV_result.txt", Results);


            //-----------------------------------------------------------------------------------------------

            //var dirString = this.rulesDirectory;
            //var dir = new DirectoryInfo(dirString);
            //dir = dir.Parent;
            //var zeoDir = new DirectoryInfo(Path.Combine(dir.FullName, "zeo++-0.3")).FullName;

            //var filepath= Path.Combine(zeoDir, "network");
            //var arguments= "-res " + Path.Combine(zeoDir, "IRMOF-1.cssr");
            //Console.WriteLine(filepath);
            //Console.WriteLine(arguments);

            //using (Process proc = new Process())
            //{
            //    // "C:\Users\zhang\AppData\Local\Programs\Julia 1.5.3\bin\\julia.exe"
            //    //C: \\Users\\zhang\\AppData\\Local\\Programs\\Julia 1.5.3\\bin\\julia.exe

            //    //C: \Users\zhang\source\repos\zeo++-0.3\network
            //    proc.StartInfo.FileName =Path.Combine(zeoDir, "network");
            //    //proc.StartInfo.Arguments = name + " -O " + name2;

            //    // "C:\Users\zhang\source\repos\MolecularSynthesis\CIFGeneration.jl"
            //    // C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\CIFGeneration.jl
            //    proc.StartInfo.Arguments = "-res " + Path.Combine(zeoDir, "IRMOF-1.cssr");
            //    //C: \Users\zhang\source\repos\MolecularSynthesis\output
            //    proc.StartInfo.WorkingDirectory = this.outputDirectory;
            //    Console.WriteLine(this.outputDirectory);
            //    //C:\\Users\\zhang\\Desktop
            //    proc.StartInfo.RedirectStandardOutput = true;
            //    proc.StartInfo.UseShellExecute=false;
            //    proc.Start();

            //    proc.WaitForExit();
            //}



            stopWatch.Stop();
            // Get the elapsed time as a TimeSpan value.
            TimeSpan ts = stopWatch.Elapsed;
            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
            ts.Hours, ts.Minutes, ts.Seconds,
            ts.Milliseconds / 10);

            Console.WriteLine("Total Running Time" + elapsedTime);

            static OBMol justMinimize(OBMol mol)
            {
                var conv = new OBConversion();
                lock (noneparallel)
                    conv.SetInAndOutFormats("pdb", "mol");

                //lock (noneparallel) 

                int ThreadNumber = System.Threading.Thread.CurrentThread.ManagedThreadId;
                string filename = "Test" + ThreadNumber.ToString() + ".mol";

                lock (noneparallel)
                    conv.WriteFile(mol, Path.Combine("/nfs/hpc/share/zhangho2/MolecularSynthesis/examples", filename));
                string minimizeOutput;

                using (Process proc = new Process())
                {
                    // C:\Program Files (x86)\OpenBabel-3.1.1

                    proc.StartInfo.FileName = "/usr/local/apps/openbabel/3.1.1/bin/obminimize";
                    proc.StartInfo.Arguments = "-ff UFF " + filename;

                    //proc.StartInfo.Arguments = "-n200 minimize.mol"; //can add arguments here like number of iterations,
                    // or '-c' convergence criteria
                    proc.StartInfo.ErrorDialog = false;
                    proc.StartInfo.WorkingDirectory = "/nfs/hpc/share/zhangho2/MolecularSynthesis/examples";
                    //proc.StartInfo.RedirectStandardError = true;
                    //proc.StartInfo.UseShellExecute = false;
                    proc.StartInfo.RedirectStandardOutput = true;
                    //proc.StartInfo.RedirectStandardInput = false;

                    Console.Write("starting OBMinimize...");
                    proc.Start();

                    minimizeOutput = proc.StandardOutput.ReadToEnd();
                    proc.WaitForExit();

                }
                lock (noneparallel)
                    conv.ReadString(mol, minimizeOutput);

                //File.Delete("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\examples\\"+ filename);




                return mol;

            }


        }
    }
}

