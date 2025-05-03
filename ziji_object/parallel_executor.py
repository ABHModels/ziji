from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import subprocess


class ParallelExecutor:
    def __init__(self, num_processes, mode, commands):
        self.num_processes = num_processes
        self.mode = mode
        self.commands = commands

    def _run_command(self, command):
        # Execute the command with no need to capture output.
        subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def execute(self):
        executor_class = (
            ProcessPoolExecutor if self.mode == "binary" else ThreadPoolExecutor
        )
        with executor_class(max_workers=self.num_processes) as executor:
            futures = [executor.submit(self._run_command, cmd) for cmd in self.commands]
            # Wait for all futures to complete (no need to process results).
            for future in futures:
                future.result()  # This will block until the command is done
