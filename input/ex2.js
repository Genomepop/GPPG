// Configuration options
{
	"individuals": 1000,
	"generations": 1000,
	"scaling": 100,
	"compression": {"name": "Greedy-Load", "k": 20, "t":10},
	//"compression": {"name": "Store-Root"},
	"genotype": {"name": "Pathway", "genes":100, "tfs":25, "regions":[10,100]},
	"operators": [{"name":"BindingSiteMutation", "overlap":5, "gainRate":1e-10, "lossRate":[1e-9,0.7] }],
	"output": {"individuals":"path/to/folder", "performance":"/Users/truths/tmp/ex2.csv"}
}