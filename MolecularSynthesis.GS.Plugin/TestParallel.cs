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
using System.Threading;
using System.Threading.Tasks;

namespace MolecularSynthesis.GS.Plugin
{
    public class TestPara : SearchProcess
    {
        // give desiredMoment
        // [] desiredMoment = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        static double[] desiredLenghtAndRadius = new double[] { 300, 50 };
        static Random rnd = new Random(0);

        public TestPara(GlobalSettings settings) : base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 2;
            AutoPlay = true;
        }


        public override string text
        {
            get { return "TestPara"; }
        }

        protected override void Run()
        {

            //TreeCandidate StartState = new TreeCandidate(seedCandidate);

            //StartState.S = 0;
            //StartState.n = 0;
            //StartState.UCB = double.MaxValue;
            //StartState.Children = new List<TreeCandidate>();

            //var option0 = rulesets[0].recognize(StartState.graph);
            //var option1 = rulesets[1].recognize(StartState.graph);
            //var option2 = rulesets[2].recognize(StartState.graph);

            ////option0[6].apply(StartState.graph, null);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[6].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[6]);

            // add a stopwatch to record time
            var timer = new Stopwatch();
            timer.Start();

            // Randomly generate .mol and .xyz files
            int TotalNumber = 100;
            var rand = new Random();
            List<string> Results = new List<string>();

            TreeCandidate StartState = new TreeCandidate(seedCandidate);

            //Parallel.For(0, number, count =>
            //{
            //    Console.WriteLine($"value of count = {count}, thread = {Thread.CurrentThread.ManagedThreadId}");
            //    //Sleep the loop for 10 miliseconds
            //    Thread.Sleep(10);
            //});


            Parallel.For(0, TotalNumber, count =>

            //for (int i = 0; i < TotalNumber; i++)
            {
                var candidate = (TreeCandidate)StartState.copy();

                var option0 = rulesets[0].recognize(candidate.graph);
                var option1 = rulesets[1].recognize(candidate.graph);
                var option2 = rulesets[2].recognize(candidate.graph);

                option0 = rulesets[0].recognize(candidate.graph);
                option0[6].apply(candidate.graph, null);
                //StartState.addToRecipe(option0[6]);

                for (int j = 0; j < 3; j++)
                {
                    //rnd.Next(0, 2); 0 or 1
                    var RuleSetNumber = rand.Next(0, 2);
                    var TotalOption = rulesets[RuleSetNumber].recognize(candidate.graph).Count;
                    var OptionNumber = rand.Next(0, TotalOption);
                    rulesets[RuleSetNumber].recognize(candidate.graph)[OptionNumber].apply(candidate.graph, null);
                    //StartState.addToRecipe(rulesets[RuleSetNumber].recognize(StartState.graph)[OptionNumber]);                    
                }

                option2 = rulesets[2].recognize(candidate.graph);
                option2[0].apply(candidate.graph, null);
                //StartState.addToRecipe(option2[0]);

                // ---------------------------------------------------------
                // code below doing parallel on HPC(about 67 threads)

                double PoreSizeValue = 0;
                // 0. recieve candidate as graph
                // 1. designgraph to mol
                var resultMol = OBFunctions.designgraphtomol(candidate.graph);
                var conv = new OBConversion();
                conv.SetInAndOutFormats("pdb", "mol");
                // generate .mol file for minimization


                int ThreadNumber = System.Threading.Thread.CurrentThread.ManagedThreadId;
                string filename1 = "Test" + ThreadNumber.ToString() + ".mol";
                Console.WriteLine("filename1:");
                Console.WriteLine(filename1);

                conv.WriteFile(resultMol, Path.Combine("/nfs/hpc/share/zhangho2/MolecularSynthesis/examples", filename1));
                string minimizeOutput;


                // 2. minimize
                using (Process proc = new Process())
                {

                    proc.StartInfo.FileName = "/usr/local/apps/openbabel/3.1.1/bin/obminimize";

                    proc.StartInfo.Arguments = "-c 1e3 -ff GAFF " + filename1;
                    //proc.StartInfo.Arguments = "-n200 minimize.mol"; //can add arguments here like number of iterations,
                    // or '-c' convergence criteria

                    //proc.StartInfo.ErrorDialog = false;
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

                conv.ReadString(resultMol, minimizeOutput);

                // after minimization, update, make new .mol file for tobacco
                OBFunctions.updatepositions(candidate.graph, resultMol);

                var FinalResultMol = OBFunctions.designgraphtomol(candidate.graph);
                conv = new OBConversion();
                conv.SetInAndOutFormats("pdb", "mol");
                // /nfs/hpc/share/zhangho2/tobacco_3.0/edges
                string filename2 = "Candidate" + ThreadNumber.ToString() + ".mol";
                Console.WriteLine("filename2:");
                Console.WriteLine(filename2);

                conv.WriteFile(FinalResultMol, Path.Combine("/nfs/hpc/share/zhangho2/MolecularSynthesis/examples", filename2));

                //File.Delete("/nfs/hpc/share/zhangho2/MolecularSynthesis/examples/TestOBD.mol");


                // 3. .mol to .cif
                
                using (Process proc = new Process())
                {

                    //C: \Users\zhang\source\repos\zeo++-0.3\network
                    // /usr/local/apps/julia/1.5/bin/julia
                    // /usr/local/apps/julia/1.5/bin/julia
                    proc.StartInfo.FileName = "/usr/local/apps/julia/1.5/bin/julia";
                    // /nfs/hpc/share/zhangho2/MolecularSynthesis


                    // make CIFGeneration.jl as function, make .mol with thread ID as input
                    proc.StartInfo.Arguments = "/nfs/hpc/share/zhangho2/MolecularSynthesis/examples/CIFGeneration.jl /nfs/hpc/share/zhangho2/MolecularSynthesis/examples/" + filename2;
                    // /nfs/hpc/share/zhangho2/tobacco_3.0/edges
                    proc.StartInfo.WorkingDirectory = "/nfs/hpc/share/zhangho2/tobacco_3.0/edges";
                    //C:\\Users\\zhang\\Desktop
                    proc.StartInfo.RedirectStandardOutput = true;
                    proc.StartInfo.UseShellExecute = false;
                    proc.Start();
                    System.Console.WriteLine("CIFGerneration.jl is running");
                    proc.WaitForExit();
                }
                Console.WriteLine("Finish .mol to .cif");

                // 3.1 move .cif to edge folder 

                var filename3 = filename2.Split(".")[0] + "_fer_tobacco.cif";
                Console.WriteLine("filename3:");
                Console.WriteLine(filename3);

                string position1 = "/nfs/hpc/share/zhangho2/MolecularSynthesis/examples/"+filename3;
                string position2 = "/nfs/hpc/share/zhangho2/tobacco_3.0/edges/" + filename3;
                System.IO.File.Move(position1, position2, true);

                // 4. invoke tobacco
                using (Process proc = new Process())
                {

                    // /usr/local/apps/python/3.8/bin/python3

                    proc.StartInfo.FileName = "/usr/local/apps/python/3.8/bin/python3";
                    // /nfs/hpc/share/zhangho2/tobacco_3.0
                    proc.StartInfo.Arguments = "/nfs/hpc/share/zhangho2/tobacco_3.0/tobacco.py";
                    //C: \Users\zhang\source\repos\MolecularSynthesis\output
                    //C: \Users\zhang\source\repos\tobacco_3.0\edges
                    proc.StartInfo.WorkingDirectory = "/nfs/hpc/share/zhangho2/tobacco_3.0";
                    proc.StartInfo.RedirectStandardOutput = true;
                    proc.StartInfo.UseShellExecute = false;
                    proc.Start();
                    System.Console.WriteLine("tabacco is running");
                    proc.WaitForExit();
                }

                // delete .cif file in edge folder to accelarate tabacco.py 
                File.Delete(position2);

                // 5. move file to the right folder to generate .cssr file
                // /nfs/hpc/share/zhangho2/tobacco_3.0/output_cifs                      

                var filename4 = "pcu_v1-6c_Zn_1_Ch_1-" + filename3;
                Console.WriteLine("filename4:");
                Console.WriteLine(filename4);

                position1 = "/nfs/hpc/share/zhangho2/tobacco_3.0/output_cifs/" + filename4;
                position2 = "/nfs/hpc/share/zhangho2/MolecularSynthesis/data/crystals/" + filename4;
                System.IO.File.Move(position1, position2, true);


                //var paths = new string[] { "/nfs", "hpc", "share", "zhangho2", "MolecularSynthesis", "data", "crystals", filename3 };
                //var NewPosition = Path.Combine(paths);
                //  /nfs/hpc/share/zhangho2/tobacco_3.0/edges
                //System.IO.File.Move("/nfs/hpc/share/zhangho2/tobacco_3.0/edges/" + filename3, NewPosition, true);

                //var folderPath = "/nfs/hpc/share/zhangho2/tobacco_3.0/output_cifs";
                //foreach (string file in Directory.EnumerateFiles(folderPath, "*.cif"))
                //{
                //    string[] details = file.Split('/');
                //    Console.WriteLine(details[7]);
                //    // /nfs/hpc/share/zhangho2/MolecularSynthesis/data/crystals
                //    var paths = new string[] { "/nfs", "hpc", "share", "zhangho2", "MolecularSynthesis", "data", "crystals", details[7] };

                //    var NewPosition = Path.Combine(paths);
                //    Console.WriteLine(NewPosition);
                //    System.IO.File.Move(file, NewPosition, true);
                //}
                Console.WriteLine("how about now???--------------------------------------------------------------");

                // invoke MakeCssr.jl

                // pcu_v1-6c_Zn_1_Ch_1-
                
                using (Process proc = new Process())
                {

                    //C: \Users\zhang\source\repos\zeo++-0.3\network
                    // /usr/bin/obminimize
                    proc.StartInfo.FileName = "/usr/local/apps/julia/1.5/bin/julia";
                    //proc.StartInfo.Arguments = name + " -O " + name2;

                    // "C:\Users\zhang\source\repos\MolecularSynthesis\CIFGeneration.jl"
                    // C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\CIFGeneration.jl
                    proc.StartInfo.Arguments = "/nfs/hpc/share/zhangho2/MolecularSynthesis/MakeCssr.jl "+ position2;
                    //C: \Users\zhang\source\repos\MolecularSynthesis\output
                    //C: \Users\zhang\source\repos\tobacco_3.0\edges
                    //C:\Users\zhang\source\repos\MolecularSynthesis\zeo++-0.3
                    proc.StartInfo.WorkingDirectory = "/nfs/hpc/share/zhangho2/MolecularSynthesis";
                    //C:\\Users\\zhang\\Desktop
                    proc.StartInfo.RedirectStandardOutput = true;
                    proc.StartInfo.UseShellExecute = false;
                    proc.Start();
                    System.Console.WriteLine("MakeCssr is running");
                    proc.WaitForExit();
                }


                // 6. invoke zeo++ to get PoresizeValue

                var filename5 = filename4.Split(".")[0] + ".cssr";
                Console.WriteLine("filename5:");
                Console.WriteLine(filename5);

                position1 = "/nfs/hpc/share/zhangho2/MolecularSynthesis/data/crystals/" + filename5;
                position2 = "/nfs/hpc/share/zhangho2/zeo++-0.3/" + filename5;
                System.IO.File.Move(position1, position2, true);


                //paths = new string[] { "/nfs", "hpc", "share", "zhangho2", "zeo++-0.3", filename5 };
                //NewPosition = Path.Combine(paths);
                //System.IO.File.Move("/nfs/hpc/share/zhangho2/zeo++-0.3/" + filename5, NewPosition, true);

                // /nfs/hpc/share/zhangho2/MolecularSynthesis
                System.Console.WriteLine("before zeo++");

                using (Process proc = new Process())
                {

                    //C: \Users\zhang\source\repos\zeo++-0.3\network
                    // /nfs/hpc/share/zhangho2/MolecularSynthesis/zeo++-0.3
                    proc.StartInfo.FileName = "/nfs/hpc/share/zhangho2/zeo++-0.3/network";
                    //proc.StartInfo.Arguments = name + " -O " + name2;

                    proc.StartInfo.Arguments = " -res " + filename5; 
                    //C: \Users\zhang\source\repos\MolecularSynthesis\output
                    proc.StartInfo.WorkingDirectory = "/nfs/hpc/share/zhangho2/zeo++-0.3";
                    //C:\\Users\\zhang\\Desktop
                    proc.StartInfo.RedirectStandardOutput = true;
                    proc.StartInfo.UseShellExecute = false;
                    proc.Start();

                    proc.WaitForExit();
                }

                //var folderPath = "/nfs/hpc/share/zhangho2/MolecularSynthesis";
                //foreach (string file in Directory.EnumerateFiles(folderPath, "*.cssr"))
                //{
                //    string[] details = file.Split('/');
                //    Console.WriteLine(details[6]);
                //    // /nfs/hpc/share/zhangho2/MolecularSynthesis/zeo++-0.3
                //    var paths = new string[] { "/nfs", "hpc", "share", "zhangho2", "zeo++-0.3", details[6] };

                //    var NewPosition = Path.Combine(paths);
                //    Console.WriteLine(NewPosition);
                //    System.IO.File.Move(file, NewPosition, true);
                //    Console.WriteLine("move .cssr file successfully");
                //    //var ZeoArgument = "-res " + details[6];
                //    //System.Console.WriteLine(file);
                //    //System.Console.WriteLine(ZeoArgument);
                //    //System.Console.WriteLine("--------------------------------------------------");

                //    using (Process proc = new Process())
                //    {

                //        //C: \Users\zhang\source\repos\zeo++-0.3\network
                //        // /nfs/hpc/share/zhangho2/MolecularSynthesis/zeo++-0.3
                //        proc.StartInfo.FileName = "/nfs/hpc/share/zhangho2/zeo++-0.3/network";
                //        //proc.StartInfo.Arguments = name + " -O " + name2;

                //        proc.StartInfo.Arguments = " -res " + details[6]; ;
                //        //C: \Users\zhang\source\repos\MolecularSynthesis\output
                //        proc.StartInfo.WorkingDirectory = "/nfs/hpc/share/zhangho2/zeo++-0.3";
                //        //C:\\Users\\zhang\\Desktop
                //        proc.StartInfo.RedirectStandardOutput = true;
                //        proc.StartInfo.UseShellExecute = false;
                //        proc.Start();

                //        proc.WaitForExit();
                //    }
                //}


                // 7. read .res file to get poresize value

                var filename6 = filename5.Split(".")[0] + ".res";
                Console.WriteLine("filename6:");
                Console.WriteLine(filename6);

                string contents = File.ReadAllText("/nfs/hpc/share/zhangho2/zeo++-0.3/" + filename6);
                string[] words = contents.Split(' ');
                Console.WriteLine(contents);

                int n = 0;
                foreach (var word in words)
                {
                    n = n + 1;
                    Console.WriteLine(n);
                    Console.WriteLine(word);
                }
                Console.WriteLine("Poresize: ", words[5]);
                //Console.WriteLine("Poresize: ", words[5]);
                //PoreSizeValue = Convert.ToDouble(words[5]);
                //Console.WriteLine("Poresizevalue: ", words[5]);

                Results.Add(PoreSizeValue.ToString());

                // /nfs/hpc/share/zhangho2/zeo++-0.3
                //folderPath = "/nfs/hpc/share/zhangho2/zeo++-0.3";
                //foreach (string file in Directory.EnumerateFiles(folderPath, "*.res"))
                //{

                //    string contents = File.ReadAllText(file);
                //    string[] words = contents.Split(' ');
                //    Console.WriteLine(contents);
                //    Console.WriteLine("Poresize: ", words[4]);

                //    //foreach (var word in words)
                //    //{
                //    //    //System.Console.WriteLine($"<{word}>");
                //    //    Console.WriteLine(word);
                //    //}
                //    PoreSizeValue = Convert.ToDouble(words[4]);
                //    Console.WriteLine("Poresizevalue: ", words[4]);
                //}

                //Results.Add(PoreSizeValue.ToString());

                //-----------------------------------------------------------
                // code below can only process one at a time

                //double PoreSize = Evaluation.GetPoreSize(candidate.graph);
                //Console.WriteLine("PoreSize= " + PoreSize);                
                //Results.Add(PoreSize.ToString());
            });

            //File.Delete("/nfs/hpc/share/zhangho2/MolecularSynthesis/examples/TestOBD.mol");


            System.IO.File.WriteAllLines(@"/nfs/hpc/share/zhangho2/MolecularSynthesis/examples_Kai/Distribution.txt", Results);


            timer.Stop();
            TimeSpan ts = timer.Elapsed;
            //string foo = "Time taken: " + timeTaken.ToString(@"m\:ss\.fff");

            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
            ts.Hours, ts.Minutes, ts.Seconds,
            ts.Milliseconds / 10);
            Console.WriteLine("RunTime " + elapsedTime);
            
        }
        
    }

}






