import logging

def configure_logging(filename='lp_solver.log', level=logging.ERROR):
    """Configures the logging system."""
    logging.basicConfig(filename=filename, level=level,
                        format='%(asctime)s - %(levelname)s - %(module)s - %(message)s')