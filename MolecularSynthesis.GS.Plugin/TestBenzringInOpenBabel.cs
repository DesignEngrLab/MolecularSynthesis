using OpenBabel;
using System;
using OpenBabelFunctions;
using GraphSynth.Search;
using GraphSynth;
using MolecularSynthesis.GS.Plugin;
using System.Xml.XPath;
using System.IO;
using System.Diagnostics;

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
            var resultMol = OBFunctions.designgraphtomol(seedGraph);
            // after convert from designgraph to .mol file
            // save the result as .mol file
            // call minimize.exe to do the energy minimization               

            resultMol = OBFunctions.InterStepMinimize(resultMol);
            OBFunctions.updatepositions(seedGraph, resultMol);
            SearchIO.addAndShowGraphWindow(seedGraph);

            var result = new double[2];
            result=Evaluation.FindLengthAndRadius(seedGraph);
            //[0] is Length, [1] is Radius, unit is 
            SearchIO.output("Length is: "+ result[0]);
            SearchIO.output("Radius is: "+ result[1]);

            var conv = new OBConversion();
            conv.SetInAndOutFormats("pdb", extension);
            conv.WriteFile(resultMol, filename + "." + extension);
            //File.AppendText()
        }
    }




}

