# Setting up and testing AmpliCI:

1.  Install Docker Desktop, which includes a command line tool that we'll use below (only do this once): [Docker Desktop](https://docs.docker.com/desktop/install/windows-install/)
1.  Go to the Windows command line.
1.  Install the AmpliCI docker image (only do this once): `docker pull kdorman/amplici:latest`
1.  Start an AmpliCI container (only do this once): `docker run -it kdorman/amplici /bin/sh`
1.  Verify the software is available, read the help: `run_AmpliCI --help`
1.  You can also try the tutorial commands as I have provided a test file (assuming you are in directory data):
    ```sh
    run_AmpliCI --fastq sim3.8.1.fastq --outfile error.out --error
    run_AmpliCI --fastq sim3.8.1.fastq --outfile test --abundance 2 --profile error.out
    run_AmpliCI --fastq sim3.8.1.fastq --outfile test.id --profile error.out --haplotypes test.fa
    ```
1.  To stop the container: `exit`
1.  Note the name of the container (last column): `docker container ls -a`


To try it out on data, the process would be generally as follows. I will type `NAME` instead of the name of the container you discovered in step 8 above.

# Running AmpliCI on your data:

1.  Go to the Windows command line
1.  Copy data to the container (assuming your data are on the `C` disk in directory `MyDataDir`): `docker cp c:\MyDataDir\MyData.fastq NAME:/data`
1.  Start the container: `docker start -a -i NAME`
1.  In the container (you should see the prompt change to `/ #`), go to where you copied the data: `cd data`
1.  Run AmpliCI commands (see the README.md), producing output file I'll call `AMPLICI_OUTPUT`.
1.  When you are done: `exit`
1.  Copy the results files to your Windows machine: `docker cp NAME:/data/AMPLICI_OUTPUT c:\MyDataDir`
