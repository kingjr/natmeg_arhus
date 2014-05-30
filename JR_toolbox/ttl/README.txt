http://tech.dir.groups.yahoo.com/group/psychtoolbox/message/11961

Re: TTL triggers in Linux
--- In psychtoolbox@yahoogroups.com, Andre Cravo <andrecravo@...> wrote:
>
> Hello,
>
> I was wondering if there is an easy way to send TTL triggers via parallel
> port in Linux (I'm currently using Ubuntu 10.10).
>
> I have seen the LabJack options, but if I understood correctly you have to
> use a USB DAQ device for it to work. I was using the " lptwrite" command in
> Windows and I was wondering if there is something similar for Linux.
>

I've just uploaded LinuxParportServer.zip to the files section of the forum. It
should work ok, although there are better solutions to come. A proper
implementation will be part of the IOPort driver in a future PTB beta.

This one works as follows:

1. unzip the file.
2. Compile the .cpp C file into a binary via (terminal window):

gcc -O2 -o parallelPortServer parallelPortServer.cpp

3. Run the executable in that terminal as root, e.g.:

sudo ./parallelPortServer

4. The server will wait for commands over a UDP network port. This is where you
use the M-File ParportTTL.m in your matlab code:

At beginnning of your script:

ParportTTL('Open', 'localhost');

At end of your script:

ParportTTL('Close', 'localhost');

For changing the state of the output lines of the parallel port in your code:

tRoundTrip = ParportTTL('Set', level [, duration=inf]);

--> See help ParportTTL for the meaning of the parameters.

This was used for a study where the stimulation computer and computer with
parallel port where different machines, therefore the network server design. You
could do the same by setting the hostname of the machine with the parallelport,
instead of 'localhost' if everything is running on the same machine.

Timing accuracy on a normal working local network is about 1 msec if i remember
correctly, but 'tRoundTrip' will give you a good upper bound on timing error for
each invocation. Last tested over 1 year ago, where the parallel port computer
was a linux box and the stimulation computer was a Apple OS/X box, so i hope it
still works "as is" on a modern system.

As i said, a proper implementation will follow as part of IOPort in the future -
one would design this quite differently nowadays, this was just quick throwaway
code for a special EEG study on a very special setup, but for the moment this
should do...

best,
-mario

> Thanks in advance
>
> --
> Andre Cravo
> Post-hoc Fellow
> University of Sao Paulo-Brazil
>



