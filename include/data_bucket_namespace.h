
/*
 Set of macros to namespace DataBucket symbols and objects.
 This is required for compatibility with PETSc v3.8
*/

#ifndef __PTATIN_DATA_BUCKET_NAMESPACE_H__
#define __PTATIN_DATA_BUCKET_NAMESPACE_H__

/* enums */
#define BTruth             ptatinBTruth
#define DataBucketViewType ptatinDataBucketViewType

/* types */
#define _p_DataField  _p_ptatinDataField
#define _p_DataBucket _p_ptatinDataBucket
#define DataField     ptatinDataField
#define DataBucket    ptatinDataBucket

/* symbols */
#define DataFieldCreate          ptatinDataFieldCreate
#define DataFieldDestroy         ptatinDataFieldDestroy
#define DataBucketCreate         ptatinDataBucketCreate
#define DataBucketDestroy        ptatinDataBucketDestroy
#define _DataBucketRegisterField _ptatinDataBucketRegisterField

#define DataFieldGetNumEntries     ptatinDataFieldGetNumEntries
#define DataFieldSetSize           ptatinDataFieldSetSize
#define DataFieldZeroBlock         ptatinDataFieldZeroBlock
#define DataFieldGetAccess         ptatinDataFieldGetAccess
#define DataFieldAccessPoint       ptatinDataFieldAccessPoint
#define DataFieldAccessPointOffset ptatinDataFieldAccessPointOffset
#define DataFieldRestoreAccess     ptatinDataFieldRestoreAccess
#define DataFieldVerifyAccess      ptatinDataFieldVerifyAccess
#define DataFieldGetAtomicSize     ptatinDataFieldGetAtomicSize

#define DataFieldGetEntries     ptatinDataFieldGetEntries
#define DataFieldRestoreEntries ptatinDataFieldRestoreEntries

#define DataFieldInsertPoint ptatinDataFieldInsertPoint
#define DataFieldCopyPoint   ptatinDataFieldCopyPoint
#define DataFieldZeroPoint   ptatinDataFieldZeroPoint

#define DataBucketGetDataFieldByName   ptatinDataBucketGetDataFieldByName
#define DataBucketQueryDataFieldByName ptatinDataBucketQueryDataFieldByName
#define DataBucketFinalize             ptatinDataBucketFinalize
#define DataBucketSetInitialSizes      ptatinDataBucketSetInitialSizes
#define DataBucketSetSizes             ptatinDataBucketSetSizes
#define DataBucketGetSizes             ptatinDataBucketGetSizes
#define DataBucketGetGlobalSizes       ptatinDataBucketGetGlobalSizes
#define DataBucketGetDataFields        ptatinDataBucketGetDataFields

#define DataBucketCopyPoint        ptatinDataBucketCopyPoint
#define DataBucketCreateFromSubset ptatinDataBucketCreateFromSubset
#define DataBucketZeroPoint        ptatinDataBucketZeroPoint

#define DataBucketLoadFromFile          ptatinDataBucketLoadFromFile
#define DataBucketView                  ptatinDataBucketView
#define DataBucketLoadRedundantFromFile ptatinDataBucketLoadRedundantFromFile

#define DataBucketAddPoint           ptatinDataBucketAddPoint
#define DataBucketRemovePoint        ptatinDataBucketRemovePoint
#define DataBucketRemovePointAtIndex ptatinDataBucketRemovePointAtIndex

#define DataBucketDuplicateFields ptatinDataBucketDuplicateFields
#define DataBucketInsertValues    ptatinDataBucketInsertValues

#define DataBucketCreatePackedArray  ptatinDataBucketCreatePackedArray
#define DataBucketDestroyPackedArray ptatinDataBucketDestroyPackedArray
#define DataBucketFillPackedArray    ptatinDataBucketFillPackedArray
#define DataBucketInsertPackedArray  ptatinDataBucketInsertPackedArray

#define DataBucketView_NATIVE          ptatinDataBucketView_NATIVE
#define DataBucketLoad_NATIVE          ptatinDataBucketLoad_NATIVE
#define DataBucketLoadRedundant_NATIVE ptatinDataBucketLoadRedundant_NATIVE

#endif
