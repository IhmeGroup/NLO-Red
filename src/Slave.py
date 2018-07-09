import Global
import logging
import Parallelism as Par
import OptimFunctions


def slave_execution():
    # MPI communicator

    MPI = Global.MPI
    comm = MPI.COMM_WORLD
    name = MPI.Get_processor_name()
    status = MPI.Status()

    rank = comm.rank
    logging.debug("Slave rank %d on %s." % (rank, name))

    # TODO put somewhere better
    tags = Par.enum('READY', 'DONE', 'EXIT', 'START', 'SLEEP', 'WAKEUP')

    # Loop until exit
    while True:
        # Notify the master the slave is ready
        comm.send(None, dest=0, tag=tags.READY)

        # Loop waiting until MPI message arrives
        while True:
            task = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
            tag = status.Get_tag()
            if not(tag == tags.WAKEUP):
                # Message received
                break
            else:
                # Already wake up
                logging.debug('Already woke up')

        # Take actions
        if tag == tags.START:
                # Start new tasks
            result = OptimFunctions.getQuantityCase(
                task[0], task[1], task[2], task[3])
            comm.send([task[4], result], dest=0, tag=tags.DONE)

        elif tag == tags.EXIT:
            # Exit signal
            break

        elif tag == tags.SLEEP:
            # Nothing, sleep until wakeup signal
            comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
            tag = status.Get_tag()
            if (tag == tags.WAKEUP):
                logging.debug('Slave waking up')
            else:
                logging.warning('There is probably something wrong')

    comm.send(None, dest=0, tag=tags.EXIT)
