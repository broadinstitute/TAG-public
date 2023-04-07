version 1.0

workflow checkBaitSetName{
  call compareBaitSetName
}

task compareBaitSetName {
  input {
    String bait_set
    File target_intervals
    File? bait_intervals
  }
  String target_intervals_name = basename(target_intervals)
  String bait_intervals_name = if defined(bait_intervals) then basename(bait_intervals) else bait_set + "."

  command <<<
python <<CODE
target_correct = ~{target_intervals_name}.startswith("~{bait_set}.")
bait_correct = ~{bait_intervals_name}.startswith("~{bait_set}.")
if not target_correct and not bait_correct:
  print("Bait and target intervals do not match the bait_set.")
  open("set_mismatch.txt", 'w').write("1")
elif not target_correct:
  print("Target intervals do not match the bait_set.")
  open("set_mismatch.txt", 'w').write("1")
elif not bait_correct:
  print("Bait intervals do not match the bait_set.")
  open("set_mismatch.txt", 'w').write("1")
else:
  print("bait_set matches the provided intervals.")
  open("set_mismatch.txt", 'w').write("0")
CODE
  >>>

  output {
    String mismatch_message = read_string(stdout())
	Int bait_mismatch = read_int("set_mismatch.txt")
  }

  runtime {
    docker: "python:3"
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    preemptible: 1
    maxRetries: 1
  }
}