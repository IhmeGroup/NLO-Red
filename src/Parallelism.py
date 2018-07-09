import Global


def enum(*sequential, **named):
    # Enumerated type in Python
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)


def tasklaunch(tasks):
    MPI = Global.MPI
    comm = MPI.COMM_WORLD   # get MPI communicator object
    size = comm.size        # total number of processes
    rank = comm.rank        # rank of this process
    status = MPI.Status()   # get MPI status object
    tags = enum('READY', 'DONE', 'EXIT', 'START', 'SLEEP', 'WAKEUP')
    partitioning = Global.partitioning

    for nt in range(comm.size-1):
        comm.send(None, dest=nt+1, tag=tags.WAKEUP)
    tasks_done = 0
    task_index = 0
    closed_workers = 0
    resultarray = [None for i in range(len(tasks))]
    launchedarray = [False for i in range(len(tasks))]

    while tasks_done < len(tasks):
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()

        if tag == tags.READY:
            foundtask = False
            indicesource = partitioning[source-1]

            for k, task in enumerate(tasks):
                if not(launchedarray[k]):
                    if task[0].isflame == False:
                        foundtask = True
                        curtask = tasks[k]
                        break
                    if task[0].partitionid in indicesource:
                        foundtask = True
                        curtask = tasks[k]
                        break

            if foundtask:
                launchedarray[k] = True
                comm.send(curtask, dest=source, tag=tags.START)
            else:
                comm.send(None, dest=source, tag=tags.SLEEP)

        elif tag == tags.DONE:
            resultarray[data[0]] = data[1]
            tasks_done += 1
            print 'Worker done', source, 'parititon', partitioning[source-1]
            print 'Task dones', tasks_done, 'Total', len(tasks)
        elif tag == tags.EXIT:
            print("Worker %d exited YOU ARE IN TROUBLE." % source)
            closed_workers += 1

    print 'Launching tasks and result collection done'
    return resultarray
