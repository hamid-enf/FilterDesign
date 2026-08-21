/*
 * fce_export_stdio.c - FILE*-backed export writer (host builds only).
 * Compile this file only when FCE_ENABLE_EXPORT_STDIO is enabled.
 */
#include "fce_internal.h"

#if FCE_ENABLE_EXPORT && FCE_ENABLE_EXPORT_STDIO

#include <stdio.h>

static size_t fce_file_writer_write(void* ctx, const char* data, size_t len)
{
    FILE* f = (FILE*)ctx;
    if (f == NULL)
        return 0;
    if (len > 0u)
        return fwrite(data, 1u, len, f);
    return 0;
}

fce_status_t fce_writer_file_init(fce_writer_t* w, void* file)
{
    if (w == NULL || file == NULL)
        return FCE_ERR_INVALID_ARGUMENT;
    w->write = fce_file_writer_write;
    w->ctx = file;
    return FCE_OK;
}

#endif /* FCE_ENABLE_EXPORT && FCE_ENABLE_EXPORT_STDIO */
