using OpenBabel;
using System;
using OpenBabelFunctions;
using GraphSynth.Search;
using GraphSynth;
using System.IO;

namespace TestOpenBabel
{
    public class TestBenzRing : SearchProcess
    {
        public override string text => "TestBenz";
        public string filename = @"..\..\..\..\test";
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
            resultMol = OBFunctions.InterStepMinimize(resultMol);
            OBFunctions.updatepositions(seedGraph, resultMol);

            var conv = new OBConversion();
            conv.SetInAndOutFormats("pdb", extension);
            conv.WriteFile(resultMol, filename + "." + extension);
            //File.AppendText()
        }
    }




}

