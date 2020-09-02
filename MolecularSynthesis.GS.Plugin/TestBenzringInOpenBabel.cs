using OpenBabel;
using System;
using OpenBabelFunctions;
using GraphSynth.Search;
using GraphSynth;
using MolecularSynthesis.Plugin;
using System.Xml.XPath;
using System.IO;
using System.Diagnostics;

namespace TestOpenBabel
{
    public class TestBenzRing : SearchProcess
    {
        public override string text => "TestBenz";
        public string filename = @"..\..\..\..\HaHa";
        public string extension = "mol";
        //deault constructor
        public TestBenzRing(GlobalSettings settings) : base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 2;
            AutoPlay = true;
        }
        protected override void Run()
        {
            var resultMol = OBFunctions.designgraphtomol(seedGraph);
            // after convert from designgraph to .mol file
            // save the result as .mol file
            // call minimize.exe to do the energy minimization
            //var conv = new OBConversion();
            //conv.SetInAndOutFormats("pdb", extension);
            //conv.WriteFile(resultMol, filename + "." + extension);
            //SearchIO.output("file has been writen");
            
            //using (Process MiniProcess = new Process())
            //{
            //    MiniProcess.StartInfo.FileName = "C:\\Program Files\\OpenBabel - 3.0.0\\obminimize.exe";
            //    MiniProcess.StartInfo.Arguments = "C:\\Users\\zhang\\source\repos\\MolecularSynthesis\\HaHa.mol";
            //    MiniProcess.Start();
                
            //}

            resultMol = OBFunctions.InterStepMinimize(resultMol);
            OBFunctions.updatepositions(seedGraph, resultMol);
            SearchIO.addAndShowGraphWindow(seedGraph);

            var result = new double[2];
            result=Evaluation.FindLengthAndRadius(seedGraph);
            //[0] is Length, [1] is Radius, unit is 
            SearchIO.output("Length is: "+ result[0]);
            SearchIO.output("Radius is: "+ result[1]);

            //var conv = new OBConversion();
            //conv.SetInAndOutFormats("pdb", extension);
            //conv.WriteFile(resultMol, filename + "." + extension);
            //File.AppendText()
        }
    }




}

