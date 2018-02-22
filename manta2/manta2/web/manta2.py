#!/usr/bin/env python2.7

import sys

import os
import sys
import string
import tempfile
import re
import glob
import shutil
#import pymongo
import time
from datetime import datetime, timedelta

from flask import Flask, request, session, g, render_template, flash, \
                  redirect, url_for, abort, send_from_directory

from werkzeug.utils import secure_filename
from flask_pymongo import PyMongo


app = Flask('manta2', static_folder='static', template_folder='templates')

#app_path = os.path.dirname(os.path.realpath(__file__))
#sys.path.append(app_path)

lib_path = os.path.abspath(os.path.join(app.root_path, os.pardir))
sys.path.append(lib_path)

import fileIO
from dbIO import search_variants

app.config.from_object('web.default_settings')
app.config.from_envvar('MANTA2_PRIVATE_SETTINGS')

app.config['ASK_FOR_FILETYPE'] = False
app.config['DID_ASK_FOR_FILETYPE'] = False
app.config['errors'] = []
app.config['warnings'] = []

mongo = PyMongo(app, config_prefix='MANTA2')


if not app.debug:
    import logging
    from logging import FileHandler

    file_handler = FileHandler(os.path.join(app.root_path, app.config['LOG_FILE']))
    file_handler.setLevel(logging.INFO)

    from logging import Formatter
    file_handler.setFormatter(Formatter(
        '%(asctime)s %(levelname)s: %(message)s '
        '[in %(pathname)s:%(lineno)d]'
    ))
    app.logger.addHandler(file_handler)


@app.teardown_appcontext
def cleanup_temp_files(error):
    """Cleanup old upload variant and download result temp dirs and files"""

    app.logger.info("Cleaning up old upload and results dirs/files")

    upload_dir = os.path.join(app.root_path, app.config['UPLOAD_REL_DIR'])
    results_dir = os.path.join(app.root_path, app.config['RESULTS_REL_DIR'])

    # Current system time in seconds
    old_file_time = datetime.now() - timedelta(days=app.config['REMOVE_RESULTFILES_OLDER_THAN'])

    for d in upload_dir, results_dir:
        for f in os.listdir(d):
            if f.startswith('tmp'):
                fpath = os.path.join(d, f)
                # File creation time in seconds
                ftime = datetime.fromtimestamp(os.path.getctime(fpath))

                if ftime < old_file_time:
                    shutil.rmtree(fpath)

    return


@app.route('/', methods=['GET'])
def start():
    return redirect(url_for('snv_input'))


@app.route('/input', methods=['GET', 'POST'])
@app.route('/upload', methods=['GET', 'POST'])
def snv_input():
    """SNV input
    """
    # Save messages locally for display and reset messages globally.
    warnings = app.config['warnings']
    app.config['warnings'] = []

    return render_template('snv_input.html',
        header='Input SNVs',
        title='Search MANTA2 for SNVs impacting TFBSs',
        ask_for_filetype=app.config['ASK_FOR_FILETYPE'],
        warnings=warnings
    )


@app.route('/results', methods=['GET', 'POST'])
def results():
    """Search the MANTA2 database for TFBSs impacted by the SNVs
    and display the results.
    """
    if request.method == 'GET':
        app.config['warnings'].append("No variant file provided.")
        return redirect(url_for('snv_input'))

    form = request.form

    app.logger.info("Got the form")

    #snv_paste = form['snv_paste']

    #snv_input = form['snv_input']
    #app.logger.info("got the snv_input")

    if 'snv_upload' not in request.files:
        err = "File upload error."
        app.logger.error(err)
        app.config['warnings'].append(err)
        return redirect(url_for('snv_input'))

    file = request.files['snv_upload']
    app.logger.info("Did the request.files")

    if file.filename == '':
        #flash('No SNV file provided')
        err = "No SNV file provided."
        app.logger.error(err)
        app.config['warnings'].append(err)
        return redirect(url_for('snv_input'))

    if not app.config['ASK_FOR_FILETYPE'] and not allowed_filename(file.filename):
        # If we didn't explicitly ask for the file type and this is not an
        # allowed extension, consider this an error. And prompt the user to
        # explicitly choose the file type.

        #flash('This is not an allowed file type.')
        msg = "This is not an allowed file type. Please provide a file with one of the following extensions ({}) or explicitly specify the file type.".format(", ".join(app.config['ALLOWED_EXTENSIONS']))

        app.logger.warn(msg)
        app.config['warnings'].append(msg)
        app.config['ASK_FOR_FILETYPE'] = True
        return redirect(url_for('snv_input'))

    snv_filename = secure_filename(file.filename)
    app.logger.info("Secure filename = {}".format(snv_filename))
    app.config['snv_filename'] = snv_filename

    upload_tempdir = tempfile.mkdtemp(dir=os.path.join(app.root_path, app.config['UPLOAD_REL_DIR'])) 

    upload_path = os.path.join(upload_tempdir, snv_filename)

    file.save(upload_path)
    #flash('File successfully uploaded')
    app.logger.info('File successfully uploaded')

    supplied_snv_filetype = None
    if app.config['ASK_FOR_FILETYPE']:
        supplied_snv_filetype = form['snv_filetype']

        # Reset back to default for next time around
        app.config['ASK_FOR_FILETYPE'] = False

        if supplied_snv_filetype:
            app.logger.info("The provided SNV filetype is {}".format(supplied_snv_filetype))
        else:
            app.config['warnings'].append("Please choose a file format.")
            return redirect(url_for('snv_input'))

    snv_filetype = None
    if supplied_snv_filetype:
        snv_filetype = supplied_snv_filetype
    else:
        # File type has not been explicitly provided. Try to determine
        # automatically.
        ft_results = fileIO.determine_file_format(upload_path)

        determined_snv_filetype = ft_results['filetype']
        ft_errors = ft_results['errors']

        if not determined_snv_filetype or ft_errors:
            # Couldn't determine the file type or there were some issues
            # determining.
            err = "Could not determine input variant file format."
            app.logger.error(err)

            app.config['warnings'].append(err)

            if ft_errors:
                for e in ft_errors:
                    app.logger.error(e)

                app.config['warnings'].extend(ft_errors)

            # Prompt for file type next time arround.
            app.config['ASK_FOR_FILETYPE'] = True

            return redirect(url_for('snv_input'))

        snv_filetype = determined_snv_filetype

    snv_filetype = snv_filetype.lower()
    app.logger.info("The SNV filetype is {}".format(snv_filetype))

    var_read_result = fileIO.read_variants_file(upload_path, snv_filetype)
    app.logger.info('Read SNVs')

    snvs = var_read_result['variants']
    var_errors = var_read_result['errors']

    if not snvs:
        no_snv_errs = []

        err = "No SNVs read from file {}.".format(snv_filename)
        app.logger.error(err)
        no_snv_errs.append(err)

        if supplied_snv_filetype:
            err = "Please make sure the file is formatted correctly according to the supplied file type of {}.".format(supplied_snv_filetype)
            app.logger.error(err)
            no_snv_errs.append(err)
        else:
            err = "Please make sure the file format matches one of the allowed file types. See the help page for a description of the file formats".format(supplied_snv_filetype)
            app.logger.error(err)
            no_snv_errs.append(err)

        if var_errors:
            # If there were some errors AND we didn't get any SNVs then
            # consider this serious and go to the error page.
            for err in var_errors:
                app.logger.error(err)

            app.config['errors'].extend(no_snv_errs)
            app.config['errors'].extend(var_errors)

            return redirect(url_for('error'))
        else:
            # Otherwise consider just some minor problems and go
            # back to the upload page.
            app.config['warnings'].extend(no_snv_errs)
            app.config['ASK_FOR_FILETYPE'] = True

            return redirect(url_for('snv_input'))

    if var_errors:
        # We had some errors but still read SNVs. In this case consider the
        # errors merely as warnings. Store then and continue.
        warning = "The following problems were encountered when reading the SNV input file:"
        app.logger.warn(warning)
        app.config['warnings'].append(warning)
        app.config['warnings'].extend(var_errors)

    results_tempdir = tempfile.mkdtemp(dir=os.path.join(app.root_path, app.config['RESULTS_REL_DIR']))

    # Full path of the results file
    results_path = os.path.join(results_tempdir, app.config['RESULTS_FILENAME'])
    #app.logger.error("results_path = {}".format(results_path))

    # Get just the temporary subdirectory name that was created above.
    temp_subdir = os.path.basename(results_tempdir)
    #app.logger.error("temp_subdir = {}".format(temp_subdir))
    
    # Relative path of the results file
    results_rel_path = os.path.join(app.config['RESULTS_REL_DIR'], temp_subdir, app.config['RESULTS_FILENAME'])
    #app.logger.error("results_rel_path = {}".format(results_rel_path))

    url_root = request.url_root
    #app.logger.error("url_root = {}".format(url_root))

    results_full_url = os.path.join(url_root, results_rel_path)
    #app.logger.error("results_full_url = {}".format(results_full_url))

    search_result = search_variants(mongo.db, snvs)

    snv_impacts = search_result['snv_impacts']
    search_errors = search_result['errors']

    if search_errors:
        # These are generally warnings more than errors.
        warning = "The following problems were encountered when searching the MANTA2 database for SNVs impacting TFBSs:"
        app.config['warnings'].append(warning)
        app.config['warnings'].extend(search_errors)


    fileIO.write_snv_impacts(results_path, snv_impacts)

    # Save warnings locally for display and clear warnings globally
    warnings = app.config['warnings']
    app.config['warnings'] = []

    return render_template(
                'results.html',
                header='Results',
                title='TFBSs impacted by SNVs',
                snv_impacts=snv_impacts,
                warnings=warnings,
                num_days=app.config['REMOVE_RESULTFILES_OLDER_THAN'],
                results_full_url=results_full_url,
                results_rel_url=results_rel_path)


@app.route('/download/results/<path:filename>', methods=['GET', 'POST'])
def download(filename):
    return send_from_directory(os.path.join(app.root_path, app.config['RESULTS_REL_DIR']), filename, as_attachment=True)


@app.route('/error', methods=['GET'])
def error():
    # Save errors locally for display and clear errors globally
    errors = app.config['errors']
    app.config['errors'] = []

    return render_template('error.html',
                            header='Error',
                            title='Error',
                            errors=errors)


#@app.route('/help', methods=['GET', 'POST'])
#def help():
#    return render_template('/help.html')


def allowed_filename(filename):
    return '.' in filename and \
        filename.rsplit('.', 1)[1].lower() in app.config['ALLOWED_EXTENSIONS']


if __name__ == '__main__':
    app.run()

